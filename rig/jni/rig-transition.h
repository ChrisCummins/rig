/*
 * Rut
 *
 * Copyright (C) 2012  Intel Corporation
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
 * License along with this library. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _RUT_TRANSITION_H_
#define _RUT_TRANSITION_H_

#include <rut.h>
#include "rig-path.h"

extern RutType rig_transition_type;

typedef struct _RigTransition RigTransition;

enum {
  RUT_TRANSITION_PROP_PROGRESS,
  RUT_TRANSITION_N_PROPS
};

/* State for an individual property that the transition is tracking */
typedef struct
{
  RutProperty *property;

  /* The transition maintains two sets of state for each property. One
   * is a constant value that is used throughout the entire transition
   * and the other is a path whose actual property value depends on
   * the progress of the timeline. Only one of these states will
   * actually be used depending on whether the property is animated.
   * However both states are retained so that if the user toggles the
   * animated button for a property information won't be lost. */

  /* path may be NULL */
  RigPath *path;
  RutBoxed constant_value;
} RigTransitionPropData;

struct _RigTransition
{
  RutObjectProps _parent;

  uint32_t id;

  float progress;

  /* Hash table of RigTransitionProperties. The key is a pointer to
   * the RutProperty indexed using g_direct_hash and the value is the
   * RigTransitionPropData struct */
  GHashTable *properties;

  RutContext *context;

  RutProperty props[RUT_TRANSITION_N_PROPS];
  RutSimpleIntrospectableProps introspectable;
};

void
rig_transition_set_progress (RigTransition *transition,
                             float progress);

RigTransitionPropData *
rig_transition_get_prop_data (RigTransition *transition,
                              RutObject *object,
                              const char *property_name);

RigPath *
rig_transition_get_path (RigTransition *transition,
                         RutObject *object,
                         const char *property_name);

RigTransition *
rig_transition_new (RutContext *context,
                    uint32_t id);

typedef void
(* RigTransitionForeachPropertyCb) (RutProperty *property,
                                    RigPath *path,
                                    const RutBoxed *constant_value,
                                    void *user_data);

void
rig_transition_foreach_property (RigTransition *transition,
                                 RigTransitionForeachPropertyCb callback,
                                 void *user_data);

void
rig_transition_free (RigTransition *transition);

#endif /* _RUT_TRANSITION_H_ */