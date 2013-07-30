/*
 * Rig
 *
 * Copyright (c) 2013 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _PARTICLE_SYSTEM_H_
#define _PARTICLE_SYSTEM_H_

#include "fuzzy.h"

/* <priv> */
struct particle_system_priv;

/*
 * A particle system
 */
struct particle_system {

	/* The type of system. */
	enum {
	  	SYSTEM_TYPE_CIRCULAR_ORBIT,
	} type;

	/* The position of the center of gravity of the system. */
	float cog[3];

	/* The standard gravitational parameter of the system center. This is the
	 * product of the gravitational constant (G) and the mass (M) of the
	 * body at the center of gravity.
	 *
	 *      μ = GM
	 */
	float u;

	/* The radius of the system. */
	struct fuzzy_float radius;

	/* The inclination of particle orbits, as an angle in radians relative
	 * to the equatorial (reference) plane. */
	struct fuzzy_float inclination;

	/* The number of particles in the system. */
	int particle_count;

	/* The size (in pixels) of particles. Each particle is represented by a
	 * rectangular point of dimensions particle_size × particle_size. */
	float particle_size;

	/* Particle color. */
	struct fuzzy_color particle_color;

	/* <priv> */
	struct particle_system_priv *priv;
};

struct particle_system *particle_system_new(CoglContext *ctx, CoglFramebuffer *fb);

void particle_system_free(struct particle_system *system);

void particle_system_paint(struct particle_system *system);

#endif /* _PARTICLE_SYSTEM_H_ */
