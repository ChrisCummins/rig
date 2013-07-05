/*
 * Rut
 *
 * Copyright (C) 2013  Intel Corporation
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

#include <config.h>

#include "rut-particle-engine.h"
#include "rut-global.h"

RutType rut_particle_engine_type;

typedef struct {
  /* Whether the particle is active or not. */
  CoglBool active;

  /* Particle velocity */
  float velocity[3];

  /* The maximum age of this particle in seconds. The particle will linearly
   * fade out until this age */
  gdouble max_age;

  /* Time to live. This value represents the age of the particle. When it
   * reaches zero, the particle ist destroyed. */
  gdouble ttl;
} ParticleEngineParticle;

typedef struct {
  float position[3];
  CoglColor color;
} ParticleEngineVertex;

struct _RutParticleEngine
{
  RutObjectProps _parent;
  int ref_count;
  RutComponentableProps component;

  RutContext *rig_ctx;

  RutMesh *mesh;
  RutMesh *pick_mesh;

  CoglBool active;

  /*
   * The maximum number of particles that can exist at any given moment in
   * time. When this number of particles has been generated, then new particles
   * will only be created as and when old particles are destroyed.
   */
  int particle_count;

  /*
   * This controls the rate at which new particles are generated.
   */
  int new_particles_per_ms;

  /*
   * The size (in pixels) of particles. Each particle is represented by a
   * rectangular point of dimensions particle_size × particle_size.
   */
  float particle_size;

  /*
   * The length of time (in seconds) that a particle exists for.
   */
  struct fuzzy_double particle_lifespan;

  /*
   * The starting position for particles.
   */
  struct fuzzy_vector particle_position;

  /*
   * A unit vector to describe particle starting direction.
   */
  struct fuzzy_vector particle_direction;

  /*
   * The initial particle speed.
   */
  struct fuzzy_float particle_speed;

  /*
   * The initial particle color. Once created, a particle maintains the same
   * color for the duration of it's lifespan, but it's opacity is related to
   * it's age, so a particle begins opaque and fades into transparency.
   */
  struct fuzzy_color particle_color;

  /*
   * A uniform global force which is applied to every particle. Can be used to
   * model gravity, wind etc.
   */
  float acceleration[3];

  /* Particle emitter priv */
  GTimer *timer;
  double current_time;
  double last_update_time;

  ParticleEngineParticle *particles;
  int active_particles_count;

  GRand *rand;

  /* Particle engine priv */
  CoglContext *ctx;
  CoglFramebuffer *fb;

  CoglPipeline *pipeline;
  CoglPrimitive *primitive;
  CoglAttributeBuffer *attribute_buffer;

  ParticleEngineVertex *vertices;

  RutSimpleIntrospectableProps introspectable;
  RutProperty properties[RUT_PARTICLE_ENGINE_N_PROPS];
};

static void
_rig_particle_engine_create_resources (RutParticleEngine *engine)
{
  CoglAttribute *attributes[2];
  unsigned int i, size = sizeof (ParticleEngineVertex);

  engine->pipeline = cogl_pipeline_new (engine->ctx);
  engine->vertices = g_new0 (ParticleEngineVertex, engine->particle_count);

  engine->attribute_buffer =
    cogl_attribute_buffer_new (engine->ctx,
                               size * engine->particle_count,
                               engine->vertices);

  attributes[0] = cogl_attribute_new (engine->attribute_buffer,
                                      "cogl_position_in", size,
                                      G_STRUCT_OFFSET (ParticleEngineVertex,
                                                       position),
                                      3, COGL_ATTRIBUTE_TYPE_FLOAT);

  attributes[1] = cogl_attribute_new (engine->attribute_buffer,
                                      "cogl_color_in", size,
                                      G_STRUCT_OFFSET (ParticleEngineVertex,
                                                       color),
                                      4, COGL_ATTRIBUTE_TYPE_FLOAT);

  engine->primitive =
    cogl_primitive_new_with_attributes (COGL_VERTICES_MODE_POINTS,
                                        engine->particle_count,
                                        attributes,
                                        G_N_ELEMENTS (attributes));

  cogl_pipeline_set_point_size (engine->pipeline, engine->particle_size);
  cogl_primitive_set_n_vertices (engine->primitive, engine->particle_count);

  for (i = 0; i < G_N_ELEMENTS (attributes); i++)
    cogl_object_unref (attributes[i]);
}

static void
_particle_engine_update(RutParticleEngine *engine)
{
  int i, updated_particles = 0, destroyed_particles = 0;
  int new_particles = 0, max_new_particles;
  gdouble tick_time;
  CoglError *error = NULL;

  /* Create resources as necessary */
  if (!engine->pipeline)
    _rig_particle_engine_create_resources (engine);

  /* Update the clocks */
  engine->last_update_time = engine->current_time;
  engine->current_time = g_timer_elapsed (engine->timer, NULL);

  tick_time = engine->current_time - engine->last_update_time;

  /* The maximum number of new particles to create for this tick. This can
   * be zero, for example in the case where the engine isn't active.
   */
  max_new_particles = engine->active ?
    tick_time * engine->new_particles_per_ms : 0;

  /* We must first map the particle engine's buffer before reading or
   * writing particle data.
   */
  engine->vertices = cogl_buffer_map (COGL_BUFFER (engine->attribute_buffer),
                                      COGL_BUFFER_ACCESS_READ_WRITE,
                                      0, &error);

  if (error != NULL) {
    g_error (G_STRLOC " failed to map buffer: %s", error->message);
    return;
  }

  /* Iterate over every particle and update/destroy/create as
   * necessary.
   */
  for (i = 0; i < engine->particle_count; i++) {
    ParticleEngineParticle *particle = &engine->particles[i];

    /* Break early if there's nothing left to do */
    if (updated_particles >= engine->active_particles_count &&
        new_particles >= max_new_particles) {
      break;
    }

    if (particle->active) {
      if (particle->ttl > 0) {
        /* Update the particle's position and color */
        //FIXME: update_particle(engine, i, tick_time);

        /* Age the particle */
        particle->ttl -= tick_time;
      } else {
        /* If a particle has expired, remove it */
        //FIXME: destroy_particle(engine, i);
        destroyed_particles++;
      }

      updated_particles++;
    } else if (new_particles < max_new_particles) {
      /* Create a particle */
      //FIXME: create_particle(engine, i);
      new_particles++;
    }
  }

  /* We can safely unmap the changes we have made to the particle buffer
   * now.
   */
  cogl_buffer_unmap(COGL_BUFFER(engine->attribute_buffer));

  /* Update particle count */
  engine->active_particles_count += new_particles - destroyed_particles;
}

static void _rut_particle_engine_free (void *object)
{
  RutParticleEngine *engine = RUT_PARTICLE_ENGINE (object);

  /* TODO: */
  g_slice_free (RutParticleEngine, engine);
}

void
_rut_particle_engine_init_type ()
{
  static RutRefCountableVTable refable_vtable = {
    rut_refable_simple_ref,
    rut_refable_simple_unref,
    _rut_particle_engine_free
  };

  static RutComponentableVTable componentable_vtable = {
    0
  };

  static RutIntrospectableVTable introspectable_vtable = {
    rut_simple_introspectable_lookup_property,
    rut_simple_introspectable_foreach_property
  };

  RutType *type = &rut_particle_engine_type;

#define TYPE RutParticleEngine

  rut_type_init (type, G_STRINGIFY (TYPE));
  rut_type_add_interface (type,
                          RUT_INTERFACE_ID_REF_COUNTABLE,
                          offsetof (TYPE, ref_count),
                          &refable_vtable);
  rut_type_add_interface (type,
                          RUT_INTERFACE_ID_COMPONENTABLE,
                          offsetof (TYPE, component),
                          &componentable_vtable);
  rut_type_add_interface (type,
                          RUT_INTERFACE_ID_INTROSPECTABLE,
                          0,
                          &introspectable_vtable);
  rut_type_add_interface (type,
                          RUT_INTERFACE_ID_SIMPLE_INTROSPECTABLE,
                          offsetof (TYPE, introspectable),
                          NULL);

#undef TYPE
}

RutParticleEngine *
rut_particle_engine_new (RutContext *ctx)
{
  RutParticleEngine *engine;

  static RutPropertySpec prop_specs[] = {
    {
      .name = "particle-engine-particle-count",
      .nick = "Particle Count",
      .type = RUT_PROPERTY_TYPE_INTEGER,
      .getter.integer_type = rut_particle_engine_get_particle_count,
      .setter.integer_type = rut_particle_engine_set_particle_count,
      .flags = RUT_PROPERTY_FLAG_READWRITE | RUT_PROPERTY_FLAG_VALIDATE,
      .validation = { .int_range = { 0, 500000 }},
      .animatable = FALSE
    },
    { NULL }
  };

  engine = g_slice_new0 (RutParticleEngine);

  rut_object_init (&engine->_parent, &rut_particle_engine_type);

  engine->ref_count = 1;
  engine->component.type = RUT_COMPONENT_TYPE_PARTICLE_ENGINE;

  engine->ctx = rut_refable_ref (ctx);

  engine->active = TRUE;

  engine->timer = g_timer_new ();
  engine->rand = g_rand_new ();

  /* TODO: defaults */

  rut_simple_introspectable_init (engine, prop_specs,
                                  engine->properties);

  return engine;
}

CoglPrimitive *
rut_particle_engine_get_primitive (RutObject *object)
{
  return NULL;
}

RutMesh *
rut_particle_engine_get_mesh (RutObject *self)
{
  return 0;
}

RutMesh *
rut_particle_engine_get_pick_mesh (RutObject *self)
{
  return 0;
}

int
rut_particle_engine_get_particle_count (RutObject *object)
{
  RutParticleEngine *engine = RUT_PARTICLE_ENGINE (object);

  return engine->particle_count;
}

void
rut_particle_engine_set_particle_count (RutObject *object,
                                        int particle_count)
{
  RutParticleEngine *engine = RUT_PARTICLE_ENGINE (object);
  const unsigned int prop = RUT_PARTICLE_ENGINE_PROP_PARTICLE_COUNT;

  if (particle_count != engine->particle_count)
    {
      engine->particle_count = particle_count;

      rut_property_dirty (&engine->rig_ctx->property_ctx,
                          &engine->properties[prop]);
    }
}
