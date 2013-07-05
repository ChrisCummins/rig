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

#ifndef __RUT_PARTICLE_ENGINE_H__
#define __RUT_PARTICLE_ENGINE_H__

#include "rut-entity.h"
#include "fuzzy.h"

typedef struct _RutParticleEngine RutParticleEngine;
#define RUT_PARTICLE_ENGINE(p) ((RutParticleEngine *)(p))
extern RutType rut_particle_engine_type;

enum {
  RUT_PARTICLE_ENGINE_PROP_PARTICLE_COUNT,

  RUT_PARTICLE_ENGINE_N_PROPS
};

void
_rut_particle_engine_init_type (void);

RutParticleEngine *
rut_particle_engine_new (RutContext *ctx);

CoglPrimitive *
rut_particle_engine_get_primitive (RutObject *object);

RutMesh *
rut_particle_engine_get_mesh (RutObject *self);

RutMesh *
rut_particle_engine_get_pick_mesh (RutObject *self);

int
rut_particle_engine_get_particle_count (RutObject *obj);

void
rut_particle_engine_set_particle_count (RutObject *obj,
                                        int particle_count);

#endif /* __RUT_PARTICLE_ENGINE_H__ */
