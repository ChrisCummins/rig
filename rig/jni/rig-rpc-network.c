/*
 * Copyright (c) 2008-2011, Dave Benson.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that the
 * following conditions are met:
 *
 * Redistributions of source code must retain the above
 * copyright notice, this list of conditions and the following
 * disclaimer.

 * Redistributions in binary form must reproduce
 * the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other
 * materials provided with the distribution.
 *
 * Neither the name
 * of "protobuf-c" nor the names of its contributors
 * may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
/*
 * Rig
 *
 * Copyright (C) 2013  Intel Corporation
 *
 * Note: Some small parts were copied from the avahi examples although
 * those files have don't name any copyright owners.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

#include <google/protobuf-c/protobuf-c-rpc.h>
#include <glib.h>

#include <rig-data.h>
#include <rig-avahi.h>

#include <rig.pb-c.h>

typedef struct _ProtobufSource
{
  GSource source;

  RigData *data;

  ProtobufCDispatch *dispatch;

  ProtobufC_RPC_Client *client;

  CoglBool pollfds_changed;

  int n_pollfds;
  GPollFD *pollfds;

  void (*connected) (RigData *data);

} ProtobufSource;

static unsigned int
protobuf_events_to_gpollfd_events (unsigned int events)
{
  return ((events & PROTOBUF_C_EVENT_READABLE) ? G_IO_IN : 0) |
    ((events & PROTOBUF_C_EVENT_WRITABLE) ? G_IO_OUT : 0);
}

static unsigned int
gpollfd_events_to_protobuf_events (unsigned int events)
{
  return ((events & G_IO_IN) ? PROTOBUF_C_EVENT_READABLE : 0) |
    ((events & G_IO_OUT) ? PROTOBUF_C_EVENT_WRITABLE : 0);
}

static int
get_timeout (ProtobufSource *protobuf_source)
{
  ProtobufCDispatch *dispatch = protobuf_source->dispatch;

  if (dispatch->has_timeout)
    {
      struct timeval tv;
      gettimeofday (&tv, NULL);

      if (dispatch->timeout_secs < (unsigned long) tv.tv_sec ||
          (dispatch->timeout_secs == (unsigned long) tv.tv_sec &&
           dispatch->timeout_usecs <= (unsigned) tv.tv_usec))
        {
          return 0;
        }
      else
        {
          int du = dispatch->timeout_usecs - tv.tv_usec;
          int ds = dispatch->timeout_secs - tv.tv_sec;
          if (du < 0)
            {
              du += 1000000;
              ds -= 1;
            }
          if (ds > INT_MAX / 1000)
            return INT_MAX / 1000 * 1000;
          else
            /* Round up, so that we ensure that something can run
               if they just wait the full duration */
            return ds * 1000 + (du + 999) / 1000;
        }
    }

  return -1;
}

static gboolean
protobuf_source_prepare (GSource *source, int *timeout)
{
  ProtobufSource *protobuf_source = (ProtobufSource *) source;
  ProtobufCDispatch *dispatch = protobuf_source->dispatch;

  if (protobuf_source->client &&
      protobuf_c_rpc_client_is_connected (protobuf_source->client))
    {
      return TRUE;
    }

  *timeout = get_timeout (protobuf_source);
  if (*timeout == 0)
    return TRUE;

  if (protobuf_source->pollfds == NULL ||
      protobuf_source->pollfds_changed ||
      dispatch->n_changes)
    {
      int i;

      /* XXX: it's possible that we could hit this path redundantly if
       * some other GSource reports that it can dispatch immediately
       * and so we will be asked to _prepare() again later and
       * dispatch->n_changes will still be set since we don't have
       * any api to clear the dispatch->changes[].
       */

      for (i = 0; i < protobuf_source->n_pollfds; i++)
        g_source_remove_poll (source, &protobuf_source->pollfds[i]);
      g_free (protobuf_source->pollfds);

      protobuf_source->n_pollfds = dispatch->n_notifies_desired;
      protobuf_source->pollfds = g_new (GPollFD, dispatch->n_notifies_desired);
      for (i = 0; i < protobuf_source->n_pollfds; i++)
        {
          protobuf_source->pollfds[i].fd = dispatch->notifies_desired[i].fd;
          protobuf_source->pollfds[i].events =
            protobuf_events_to_gpollfd_events (dispatch->notifies_desired[i].events);

          g_source_add_poll (source, &protobuf_source->pollfds[i]);
        }
    }

  protobuf_source->pollfds_changed = FALSE;

  return FALSE;
}

static gboolean
protobuf_source_check (GSource *source)
{
  ProtobufSource *protobuf_source = (ProtobufSource *) source;
  ProtobufCDispatch *dispatch = protobuf_source->dispatch;
  int i;

  /* XXX: when we call protobuf_c_dispatch_dispatch() that will clear
   * dispatch->changes[] and so we make sure to check first if there
   * have been changes made to the pollfds so later when we prepare
   * to poll again we can update the GSource pollfds.
   */
  if (dispatch->n_changes)
    protobuf_source->pollfds_changed = TRUE;

  if (dispatch->has_idle)
    return TRUE;

  for (i = 0; i < protobuf_source->n_pollfds; i++)
    if (protobuf_source->pollfds[i].revents)
      return TRUE;

  if (get_timeout (protobuf_source) == 0)
    return TRUE;

  return FALSE;
}

static gboolean
protobuf_source_dispatch (GSource *source,
                          GSourceFunc callback,
                          void *user_data)
{
  ProtobufSource *protobuf_source = (ProtobufSource *) source;
  GPollFD *gpollfds = protobuf_source->pollfds;
  int i;
  int n_events;
  ProtobufC_FDNotify *events;
  void *to_free = NULL;


  if (protobuf_source->client &&
      protobuf_c_rpc_client_is_connected (protobuf_source->client))
    {
      if (protobuf_source->connected)
        {
          protobuf_source->connected (protobuf_source->data);
          protobuf_source->connected = NULL;
        }
    }

  n_events = 0;
  for (i = 0; i < protobuf_source->n_pollfds; i++)
    if (gpollfds[i].revents)
      n_events++;

  if (n_events < 128)
    events = alloca (sizeof (ProtobufC_FDNotify) * n_events);
  else
    to_free = events = g_new (ProtobufC_FDNotify, n_events);

  n_events = 0;
  for (i = 0; i < protobuf_source->n_pollfds; i++)
    if (gpollfds[i].revents)
      {
        events[n_events].fd = gpollfds[i].fd;
        events[n_events].events =
          gpollfd_events_to_protobuf_events (gpollfds[i].revents);

        /* note that we may actually wind up with fewer events after
         * calling gpollfd_events_to_protobuf_events() */
        if (events[n_events].events != 0)
          n_events++;
      }

  protobuf_c_dispatch_dispatch (protobuf_source->dispatch, n_events, events);

  if (to_free)
    g_free (to_free);

  return TRUE;
}

static GSourceFuncs
protobuf_source_funcs =
  {
    protobuf_source_prepare,
    protobuf_source_check,
    protobuf_source_dispatch,
    NULL
  };

static GSource *
protobuf_source_new (RigData *data)
{
  GSource *source = g_source_new (&protobuf_source_funcs,
                                  sizeof (ProtobufSource));
  ProtobufSource *protobuf_source = (ProtobufSource *)source;

  protobuf_source->dispatch = protobuf_c_dispatch_default ();
  protobuf_source->data = data;

  return source;
}

static void
example__test (Rig__Slave_Service *service,
               const Rig__Query *query,
               Rig__TestResult_Closure closure,
               void *closure_data)
{
  Rig__TestResult result = RIG__TEST_RESULT__INIT;

  g_return_if_fail (query != NULL);

  closure (&result, closure_data);
}

static Rig__Slave_Service rig_slave_service =
  RIG__SLAVE__INIT(example__);

void
rig_start_rpc_server (RigData *data, int port)
{
  char *port_str = g_strdup_printf ("%d", port);
  GSource *source;

  data->rpc_server =
    protobuf_c_rpc_server_new (PROTOBUF_C_RPC_ADDRESS_TCP,
                               port_str,
                               (ProtobufCService *) &rig_slave_service,
                               NULL);

  g_free (port_str);

  source = protobuf_source_new (data);
  g_source_attach (source, NULL);

  data->network_port = port;

  rig_avahi_register_service (data);
}

void
rig_start_rpc_client (RigData *data,
                      const char *hostname,
                      int port,
                      void (*connected) (RigData *data))
{
  char *addr_str = g_strdup_printf ("%s:%d", hostname, port);
  ProtobufC_RPC_Client *client;
  GSource *source;
  ProtobufSource *protobuf_source;

  client =
    protobuf_c_rpc_client_new (PROTOBUF_C_RPC_ADDRESS_TCP,
                               addr_str,
                               &rig__slave__descriptor,
                               NULL);

  g_free (addr_str);

  source = protobuf_source_new (data);
  protobuf_source = (ProtobufSource *)source;

  protobuf_source->client = client;
  protobuf_source->connected = connected;

  g_source_attach (source, NULL);

  data->network_port = port;

  data->rpc_client = client;
}