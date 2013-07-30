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

#ifndef _PARTICLE_ENGINE_H_
#define _PARTICLE_ENGINE_H_

#include <cogl/cogl.h>

/*
 * The particle engine is an opaque data structure
 */
struct particle_engine;

/*
 * Create and return a new particle engine.
 */
struct particle_engine *particle_engine_new(CoglContext *ctx,
					    CoglFramebuffer *fb,
					    int particle_count,
					    float particle_size);

/*
 * Destroy a particle engine and free associated resources.
 */
void particle_engine_free(struct particle_engine *engine);

/*
 * This maps the particle's vertices buffer according to the given access
 * flags.
 */
inline void particle_engine_push_buffer(struct particle_engine *engine,
					CoglBufferAccess access,
					CoglBufferMapHint hints);

/*
 * This unmaps the internal attribute buffer, writing out any changes made to
 * the particle vertices.
 */
inline void particle_engine_pop_buffer(struct particle_engine *engine);

/*
 * Returns a pointer to the given particle's position as an array of floats [x, y, z].
 */
inline float *particle_engine_get_particle_position(struct particle_engine *engine, int index);

/*
 * Returns a pointer to the given particle's color.
 */
inline CoglColor *particle_engine_get_particle_color(struct particle_engine *engine, int index);

/*
 * Paint function.
 */
void particle_engine_paint(struct particle_engine *engine);

#endif /* _PARTICLE_ENGINE_H */
