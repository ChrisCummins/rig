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
 * License along with this library. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <glib.h>

#include <cogl/cogl.h>

#include "rut-camera.h"
#include "rut-context.h"
#include "rut-global.h"
#include "rut-camera-private.h"

typedef struct
{
  float x, y, z, w;
  float r, g, b, a;
} RutVertex4C4;

static CoglUserDataKey fb_camera_key;

static RutPropertySpec _rut_camera_prop_specs[] = {
  {
    .name = "mode",
    .nick = "Mode",
    .type = RUT_PROPERTY_TYPE_ENUM,
    .getter = rut_camera_get_projection_mode,
    .setter = rut_camera_set_projection_mode,
    .flags = RUT_PROPERTY_FLAG_READWRITE |
      RUT_PROPERTY_FLAG_VALIDATE,
    .validation = { .ui_enum = &_rut_projection_ui_enum }
  },
  {
    .name = "viewport_x",
    .nick = "Viewport X",
    .type = RUT_PROPERTY_TYPE_FLOAT,
    .data_offset = offsetof (RutCamera, viewport[0]),
    .setter = rut_camera_set_viewport_x
  },
  {
    .name = "viewport_y",
    .nick = "Viewport Y",
    .type = RUT_PROPERTY_TYPE_FLOAT,
    .data_offset = offsetof (RutCamera, viewport[1]),
    .setter = rut_camera_set_viewport_y
  },
  {
    .name = "viewport_width",
    .nick = "Viewport Width",
    .type = RUT_PROPERTY_TYPE_FLOAT,
    .data_offset = offsetof (RutCamera, viewport[2]),
    .setter = rut_camera_set_viewport_width
  },
  {
    .name = "viewport_height",
    .nick = "Viewport Height",
    .type = RUT_PROPERTY_TYPE_FLOAT,
    .data_offset = offsetof (RutCamera, viewport[3]),
    .setter = rut_camera_set_viewport_height
  },
  {
    .name = "fov",
    .nick = "Field Of View",
    .type = RUT_PROPERTY_TYPE_INTEGER,
    .getter = rut_camera_get_field_of_view,
    .setter = rut_camera_set_field_of_view,
    .flags = RUT_PROPERTY_FLAG_READWRITE |
      RUT_PROPERTY_FLAG_VALIDATE,
    .validation = { .int_range = { .min = 1, .max = 135}},
    .animatable = TRUE
  },
  {
    .name = "near",
    .nick = "Near Plane",
    .type = RUT_PROPERTY_TYPE_FLOAT,
    .getter = rut_camera_get_near_plane,
    .setter = rut_camera_set_near_plane,
    .flags = RUT_PROPERTY_FLAG_READWRITE,
    .animatable = TRUE
  },
  {
    .name = "far",
    .nick = "Far Plane",
    .type = RUT_PROPERTY_TYPE_FLOAT,
    .getter = rut_camera_get_far_plane,
    .setter = rut_camera_set_far_plane,
    .flags = RUT_PROPERTY_FLAG_READWRITE,
    .animatable = TRUE
  },
  {
    .name = "background_color",
    .nick = "Background Color",
    .type = RUT_PROPERTY_TYPE_COLOR,
    .getter = rut_camera_get_background_color,
    .setter = rut_camera_set_background_color,
    .flags = RUT_PROPERTY_FLAG_READWRITE,
    .animatable = TRUE
  },
  /* FIXME: Figure out how to expose the orthographic coordinates as
   * properties? */
  { 0 }
};


static void
_rut_camera_free (void *object)
{
  RutCamera *camera = object;
  GList *l;

  cogl_object_unref (camera->fb);

  for (l = camera->input_regions; l; l = l->next)
    rut_refable_unref (l->data);
  g_list_free (camera->input_regions);

  g_slice_free (RutCamera, object);
}

RutRefCountableVTable _rut_camera_ref_countable = {
  rut_refable_simple_ref,
  rut_refable_simple_unref,
  _rut_camera_free
};

#if 0
static RutGraphableVTable _rut_camera_graphable_vtable = {
  NULL, /* child remove */
  NULL, /* child add */
  NULL /* parent changed */
};

static RutTransformableVTable _rut_camera_transformable_vtable = {
  rut_camera_get_view_transform
};
#endif

CoglPrimitive *
rut_camera_create_frustum_primitive (RutCamera *camera)
{
  RutVertex4C4 vertices[8] = {
    /* near plane in projection space */
    {-1, -1, -1, 1, /* position */ .8, .8, .8, .8 /* color */ },
    { 1, -1, -1, 1,                .8, .8, .8, .8             },
    { 1,  1, -1, 1,                .8, .8, .8, .8             },
    {-1,  1, -1, 1,                .8, .8, .8, .8             },
    /* far plane in projection space */
    {-1, -1, 1, 1,  /* position */ .3, .3, .3, .3 /* color */ },
    { 1, -1, 1, 1,                 .3, .3, .3, .3             },
    { 1,  1, 1, 1,                 .3, .3, .3, .3             },
    {-1,  1, 1, 1,                 .3, .3, .3, .3             }
  };
  const CoglMatrix *projection_inv;
  CoglAttributeBuffer *attribute_buffer;
  CoglAttribute *attributes[2];
  CoglPrimitive *primitive;
  CoglIndices *indices;
  unsigned char indices_data[24] = {
      0,1, 1,2, 2,3, 3,0,
      4,5, 5,6, 6,7, 7,4,
      0,4, 1,5, 2,6, 3,7
  };
  int i;

  projection_inv = rut_camera_get_inverse_projection (camera);

  for (i = 0; i < 8; i++)
    {
      cogl_matrix_transform_point (projection_inv,
                                   &vertices[i].x,
                                   &vertices[i].y,
                                   &vertices[i].z,
                                   &vertices[i].w);
      vertices[i].x /= vertices[i].w;
      vertices[i].y /= vertices[i].w;
      vertices[i].z /= vertices[i].w;
      vertices[i].w /= 1.0f;
    }

  attribute_buffer = cogl_attribute_buffer_new (rut_cogl_context,
                                                8 * sizeof (RutVertex4C4),
                                                vertices);

  attributes[0] = cogl_attribute_new (attribute_buffer,
                                      "cogl_position_in",
                                      sizeof (RutVertex4C4),
                                      offsetof (RutVertex4C4, x),
                                      3,
                                      COGL_ATTRIBUTE_TYPE_FLOAT);

  attributes[1] = cogl_attribute_new (attribute_buffer,
                                      "cogl_color_in",
                                      sizeof (RutVertex4C4),
                                      offsetof (RutVertex4C4, r),
                                      4,
                                      COGL_ATTRIBUTE_TYPE_FLOAT);

  indices = cogl_indices_new (rut_cogl_context,
                              COGL_INDICES_TYPE_UNSIGNED_BYTE,
                              indices_data,
                              G_N_ELEMENTS (indices_data));

  primitive = cogl_primitive_new_with_attributes (COGL_VERTICES_MODE_LINES,
                                                  8, attributes, 2);

  cogl_primitive_set_indices (primitive, indices, G_N_ELEMENTS(indices_data));

  cogl_object_unref (attribute_buffer);
  cogl_object_unref (attributes[0]);
  cogl_object_unref (attributes[1]);
  cogl_object_unref (indices);

  return primitive;
}

typedef struct _CameraFlushState
{
  RutCamera *current_camera;
  unsigned int transform_age;
  CoglBool in_frame;
} CameraFlushState;

static void
free_camera_flush_state (void *user_data)
{
  CameraFlushState *state = user_data;
  g_slice_free (CameraFlushState, state);
}

static void
_rut_camera_flush_transforms (RutCamera *camera)
{
  const CoglMatrix *projection;
  CoglFramebuffer *fb = camera->fb;
  CameraFlushState *state =
    cogl_object_get_user_data (COGL_OBJECT (fb), &fb_camera_key);
  if (!state)
    {
      state = g_slice_new (CameraFlushState);
      state->in_frame = FALSE;
      cogl_object_set_user_data (COGL_OBJECT (fb),
                                 &fb_camera_key,
                                 state,
                                 free_camera_flush_state);
    }
  else if (state->current_camera == camera &&
           camera->transform_age == state->transform_age)
    {
      state->in_frame = TRUE;
      return;
    }

  if (state->in_frame)
    {
      g_warning ("Un-balanced rut_camera_flush/_end calls: "
                 "repeat _flush() calls before _end()");
    }

  cogl_framebuffer_set_viewport (fb,
                                 camera->viewport[0],
                                 camera->viewport[1],
                                 camera->viewport[2],
                                 camera->viewport[3]);

  projection = rut_camera_get_projection (camera);
  cogl_framebuffer_set_projection_matrix (fb, projection);

  cogl_framebuffer_set_modelview_matrix (fb, &camera->view);

  state->current_camera = camera;
  state->transform_age = camera->transform_age;
  state->in_frame = TRUE;
}

static RutComponentableVTable _rut_camera_componentable_vtable = {
  .draw = NULL
};

static RutIntrospectableVTable _rut_camera_introspectable_vtable = {
  rut_simple_introspectable_lookup_property,
  rut_simple_introspectable_foreach_property
};


RutType rut_camera_type;

void
_rut_camera_init_type (void)
{
  rut_type_init (&rut_camera_type);
  rut_type_add_interface (&rut_camera_type,
                           RUT_INTERFACE_ID_REF_COUNTABLE,
                           offsetof (RutCamera, ref_count),
                           &_rut_camera_ref_countable);
  rut_type_add_interface (&rut_camera_type,
                           RUT_INTERFACE_ID_COMPONENTABLE,
                           offsetof (RutCamera, component),
                           &_rut_camera_componentable_vtable);
  rut_type_add_interface (&rut_camera_type,
                          RUT_INTERFACE_ID_INTROSPECTABLE,
                          0, /* no implied properties */
                          &_rut_camera_introspectable_vtable);
  rut_type_add_interface (&rut_camera_type,
                          RUT_INTERFACE_ID_SIMPLE_INTROSPECTABLE,
                          offsetof (RutCamera, introspectable),
                          NULL); /* no implied vtable */

#if 0
  rut_type_add_interface (&rut_camera_type,
                           RUT_INTERFACE_ID_GRAPHABLE,
                           offsetof (RutCamera, graphable),
                           &_rut_camera_graphable_vtable);
  rut_type_add_interface (&rut_camera_type,
                           RUT_INTERFACE_ID_TRANSFORMABLE,
                           0,
                           &_rut_camera_transformable_vtable);
#endif
}

RutCamera *
rut_camera_new (RutContext *ctx, CoglFramebuffer *framebuffer)
{
  RutCamera *camera = g_slice_new0 (RutCamera);
  int width = cogl_framebuffer_get_width (framebuffer);
  int height = cogl_framebuffer_get_height (framebuffer);

  camera->ctx = rut_refable_ref (ctx);

  rut_object_init (&camera->_parent, &rut_camera_type);

  camera->ref_count = 1;

  camera->component.type = RUT_COMPONENT_TYPE_CAMERA;

  rut_camera_set_background_color4f (camera, 0, 0, 0, 1);
  camera->clear_fb = TRUE;

  //rut_graphable_init (RUT_OBJECT (camera));

  camera->orthographic = TRUE;
  camera->x1 = 0;
  camera->y1 = 0;
  camera->x2 = width;
  camera->y2 = height;
  camera->near = -1;
  camera->far = 100;

  camera->projection_cache_age = -1;
  camera->inverse_projection_age = -1;

  cogl_matrix_init_identity (&camera->view);
  camera->inverse_view_age = -1;

  camera->viewport[2] = width;
  camera->viewport[3] = height;

  camera->transform_age = 0;

  cogl_matrix_init_identity (&camera->input_transform);

  camera->fb = cogl_object_ref (framebuffer);

  camera->frame = 0;
  camera->timer = g_timer_new ();
  g_timer_start (camera->timer);

  rut_simple_introspectable_init (camera,
                                  _rut_camera_prop_specs,
                                  camera->properties);

  return camera;
}

void
rut_camera_set_background_color4f (RutCamera *camera,
                                   float red,
                                   float green,
                                   float blue,
                                   float alpha)
{
  cogl_color_init_from_4f (&camera->bg_color,
                           red, green, blue, alpha);
  rut_property_dirty (&camera->ctx->property_ctx,
                      &camera->properties[RUT_CAMERA_PROP_BG_COLOR]);
}

void
rut_camera_set_background_color (RutCamera *camera,
                                 const CoglColor *color)
{
  camera->bg_color = *color;
  rut_property_dirty (&camera->ctx->property_ctx,
                      &camera->properties[RUT_CAMERA_PROP_BG_COLOR]);
}

void
rut_camera_get_background_color (RutCamera *camera,
                                 CoglColor *color)
{
  *color = camera->bg_color;
}

void
rut_camera_set_clear (RutCamera *camera,
                      CoglBool clear)
{
  if (clear)
    camera->clear_fb = TRUE;
  else
    camera->clear_fb = FALSE;
}

CoglFramebuffer *
rut_camera_get_framebuffer (RutCamera *camera)
{
  return camera->fb;
}

void
rut_camera_set_framebuffer (RutCamera *camera,
                            CoglFramebuffer *framebuffer)
{
  if (camera->fb == framebuffer)
    return;

  cogl_object_unref (camera->fb);
  camera->fb = cogl_object_ref (framebuffer);
}

static void
_rut_camera_set_viewport (RutCamera *camera,
                          float x,
                          float y,
                          float width,
                          float height)
{
  if (camera->viewport[0] == x &&
      camera->viewport[1] == y &&
      camera->viewport[2] == width &&
      camera->viewport[3] == height)
    return;

  /* If the aspect ratio changes we may need to update the projection
   * matrix... */
  if ((!camera->orthographic) &&
      (camera->viewport[2] / camera->viewport[3]) != (width / height))
    camera->projection_age++;

  camera->viewport[0] = x;
  camera->viewport[1] = y;
  camera->viewport[2] = width;
  camera->viewport[3] = height;

  camera->transform_age++;
}

void
rut_camera_set_viewport (RutCamera *camera,
                         float x,
                         float y,
                         float width,
                         float height)
{
  _rut_camera_set_viewport (camera, x, y, width, height);
  rut_property_dirty (&camera->ctx->property_ctx,
                      &camera->properties[RUT_CAMERA_PROP_VIEWPORT_X]);
  rut_property_dirty (&camera->ctx->property_ctx,
                      &camera->properties[RUT_CAMERA_PROP_VIEWPORT_Y]);
  rut_property_dirty (&camera->ctx->property_ctx,
                      &camera->properties[RUT_CAMERA_PROP_VIEWPORT_WIDTH]);
  rut_property_dirty (&camera->ctx->property_ctx,
                      &camera->properties[RUT_CAMERA_PROP_VIEWPORT_HEIGHT]);
}

void
rut_camera_set_viewport_x (RutCamera *camera,
                           float x)
{
  _rut_camera_set_viewport (camera,
                            x,
                            camera->viewport[1],
                            camera->viewport[2],
                            camera->viewport[3]);
  rut_property_dirty (&camera->ctx->property_ctx,
                      &camera->properties[RUT_CAMERA_PROP_VIEWPORT_X]);
}

void
rut_camera_set_viewport_y (RutCamera *camera,
                           float y)
{
  _rut_camera_set_viewport (camera,
                            camera->viewport[0],
                            y,
                            camera->viewport[2],
                            camera->viewport[3]);
  rut_property_dirty (&camera->ctx->property_ctx,
                      &camera->properties[RUT_CAMERA_PROP_VIEWPORT_Y]);
}

void
rut_camera_set_viewport_width (RutCamera *camera,
                               float width)
{
  _rut_camera_set_viewport (camera,
                            camera->viewport[0],
                            camera->viewport[1],
                            width,
                            camera->viewport[3]);
  rut_property_dirty (&camera->ctx->property_ctx,
                      &camera->properties[RUT_CAMERA_PROP_VIEWPORT_WIDTH]);
}

void
rut_camera_set_viewport_height (RutCamera *camera,
                                float height)
{
  _rut_camera_set_viewport (camera,
                            camera->viewport[0],
                            camera->viewport[1],
                            camera->viewport[2],
                            height);
  rut_property_dirty (&camera->ctx->property_ctx,
                      &camera->properties[RUT_CAMERA_PROP_VIEWPORT_HEIGHT]);
}

const float *
rut_camera_get_viewport (RutCamera *camera)
{
  return camera->viewport;
}

const CoglMatrix *
rut_camera_get_projection (RutCamera *camera)
{
  if (G_UNLIKELY (camera->projection_cache_age != camera->projection_age))
    {
      cogl_matrix_init_identity (&camera->projection);

      if (camera->orthographic)
        {
          cogl_matrix_orthographic (&camera->projection,
                                    camera->x1,
                                    camera->y1,
                                    camera->x2,
                                    camera->y2,
                                    camera->near,
                                    camera->far);
        }
      else
        {
          float aspect_ratio = camera->viewport[2] / camera->viewport[3];
          cogl_matrix_perspective (&camera->projection,
                                   camera->fov,
                                   aspect_ratio,
                                   camera->near,
                                   camera->far);
        }

      camera->projection_cache_age = camera->projection_age;
    }

  return &camera->projection;
}

void
rut_camera_set_near_plane (RutCamera *camera,
                           float near)
{
  if (camera->near == near)
    return;

  camera->near = near;
  rut_property_dirty (&camera->ctx->property_ctx,
                      &camera->properties[RUT_CAMERA_PROP_NEAR]);
  camera->projection_age++;
  camera->transform_age++;
}

float
rut_camera_get_near_plane (RutCamera *camera)
{
  return camera->near;
}

void
rut_camera_set_far_plane (RutCamera *camera,
                          float far)
{
  if (camera->far == far)
    return;

  camera->far = far;
  rut_property_dirty (&camera->ctx->property_ctx,
                      &camera->properties[RUT_CAMERA_PROP_FAR]);
  camera->projection_age++;
  camera->transform_age++;
}

float
rut_camera_get_far_plane (RutCamera *camera)
{
  return camera->far;
}

RutProjection
rut_camera_get_projection_mode (RutCamera *camera)
{
  if (camera->orthographic)
    return RUT_PROJECTION_ORTHOGRAPHIC;
  else
    return RUT_PROJECTION_PERSPECTIVE;
}

void
rut_camera_set_projection_mode (RutCamera  *camera,
                                RutProjection projection)
{
  CoglBool orthographic;

  if (projection == RUT_PROJECTION_ORTHOGRAPHIC)
    orthographic = TRUE;
  else
    orthographic = FALSE;

  if (orthographic != camera->orthographic)
    {
      camera->orthographic = orthographic;
      rut_property_dirty (&camera->ctx->property_ctx,
                          &camera->properties[RUT_CAMERA_PROP_MODE]);
      camera->projection_age++;
      camera->transform_age++;
    }
}

void
rut_camera_set_field_of_view (RutCamera *camera,
                              float fov)
{
  if (camera->fov == fov)
    return;

  camera->fov = fov;
  rut_property_dirty (&camera->ctx->property_ctx,
                      &camera->properties[RUT_CAMERA_PROP_FOV]);
  if (!camera->orthographic)
    {
      camera->projection_age++;
      camera->transform_age++;
    }
}

float
rut_camera_get_field_of_view (RutCamera *camera)
{
  return camera->fov;
}

void
rut_camera_set_orthographic_coordinates (RutCamera *camera,
                                         float x1,
                                         float y1,
                                         float x2,
                                         float y2)
{
  if (camera->x1 == x1 &&
      camera->y1 == y1 &&
      camera->x2 == x2 &&
      camera->y2 == y2)
    return;

  camera->x1 = x1;
  camera->y1 = y1;
  camera->x2 = x2;
  camera->y2 = y2;

  if (camera->orthographic)
    camera->projection_age++;
}

const CoglMatrix *
rut_camera_get_inverse_projection (RutCamera *camera)
{
  const CoglMatrix *projection;

  if (camera->inverse_projection_age == camera->projection_age)
    return &camera->inverse_projection;

  projection = rut_camera_get_projection (camera);

  if (!cogl_matrix_get_inverse (projection,
                                &camera->inverse_projection))
    return NULL;

  camera->inverse_projection_age = camera->projection_age;
  return &camera->inverse_projection;
}

void
rut_camera_set_view_transform (RutCamera *camera,
                               const CoglMatrix *view)
{
  camera->view = *view;

  camera->view_age++;
  camera->transform_age++;

  /* XXX: we have no way to assert that we are at the bottom of the
   * matrix stack at this point, so this might do bad things...
   */
  //cogl_framebuffer_set_modelview_matrix (camera->fb,
  //                                       &camera->view);
}

const CoglMatrix *
rut_camera_get_view_transform (RutCamera *camera)
{
  return &camera->view;
}

const CoglMatrix *
rut_camera_get_inverse_view_transform (RutCamera *camera)
{
  if (camera->inverse_view_age == camera->view_age)
    return &camera->inverse_view;

  if (!cogl_matrix_get_inverse (&camera->view,
                                &camera->inverse_view))
    return NULL;

  camera->inverse_view_age = camera->view_age;
  return &camera->inverse_view;
}

void
rut_camera_set_input_transform (RutCamera *camera,
                                const CoglMatrix *input_transform)
{
  camera->input_transform = *input_transform;
}

void
rut_camera_add_input_region (RutCamera *camera,
                             RutInputRegion *region)
{
  g_print ("add input region %p, %p\n", camera, region);
  camera->input_regions = g_list_prepend (camera->input_regions, region);
}

void
rut_camera_remove_input_region (RutCamera *camera,
                                RutInputRegion *region)
{
  camera->input_regions = g_list_remove (camera->input_regions, region);
}

#if 0
typedef struct _RutCameraInputCallbackState
{
  RutInputCallback callback;
  void *user_data;
} RutCameraInputCallbackState;

RutInputEventStatus
_rut_camera_input_callback_wrapper (RutCameraInputCallbackState *state,
                                    RutInputEvent *event)
{
  RutCamera *camera = state->camera;
  float *viewport = camera->viewport;

  cogl_matrix_translate (event->input_transform,
                         viewport[0], viewport[1], 0);
  cogl_matrix_scale (event->input_transform,
                     cogl_framebuffer_get_width (camera->fb) / viewport[2],
                     cogl_framebuffer_get_height (camera->fb) / viewport[3]);

  return state->callback (event, state->user_data);
}

void
rut_camera_add_input_callback (RutCamera *camera,
                               RutInputCallback callback,
                               void *user_data)
{
  RutCameraInputCallbackState *state = g_slice_new (RutCameraInputCallbackState);
  state->camera = camera;
  state->callback = callback;
  state->user_data = user_data;
  camera->input_callbacks = g_list_prepend (camera->input_callbacks, state);
}
#endif

CoglBool
rut_camera_transform_window_coordinate (RutCamera *camera,
                                        float *x,
                                        float *y)
{
  float *viewport = camera->viewport;
  *x -= viewport[0];
  *y -= viewport[1];

  if (*x < 0 || *x >= viewport[2] || *y < 0 || *y >= viewport[3])
    return FALSE;
  else
    return TRUE;
}

void
rut_camera_unproject_coord (RutCamera *camera,
                            const CoglMatrix *modelview,
                            const CoglMatrix *inverse_modelview,
                            float object_coord_z,
                            float *x,
                            float *y)
{
  const CoglMatrix *projection = rut_camera_get_projection (camera);
  const CoglMatrix *inverse_projection =
    rut_camera_get_inverse_projection (camera);
  //float z;
  float ndc_x, ndc_y, ndc_z, ndc_w;
  float eye_x, eye_y, eye_z, eye_w;
  const float *viewport = rut_camera_get_viewport (camera);

  /* Convert item z into NDC z */
  {
    //float x = 0, y = 0, z = 0, w = 1;
    float z = 0, w = 1;
    float tmp_x, tmp_y, tmp_z;
    const CoglMatrix *m = modelview;

    tmp_x = m->xw;
    tmp_y = m->yw;
    tmp_z = m->zw;

    m = projection;
    z = m->zx * tmp_x + m->zy * tmp_y + m->zz * tmp_z + m->zw;
    w = m->wx * tmp_x + m->wy * tmp_y + m->wz * tmp_z + m->ww;

    ndc_z = z / w;
  }

  /* Undo the Viewport transform, putting us in Normalized Device Coords */
  ndc_x = (*x - viewport[0]) * 2.0f / viewport[2] - 1.0f;
  ndc_y = ((viewport[3] - 1 + viewport[1] - *y) * 2.0f / viewport[3] - 1.0f);

  /* Undo the Projection, putting us in Eye Coords. */
  ndc_w = 1;
  cogl_matrix_transform_point (inverse_projection,
                               &ndc_x, &ndc_y, &ndc_z, &ndc_w);
  eye_x = ndc_x / ndc_w;
  eye_y = ndc_y / ndc_w;
  eye_z = ndc_z / ndc_w;
  eye_w = 1;

  /* Undo the Modelview transform, putting us in Object Coords */
  cogl_matrix_transform_point (inverse_modelview,
                               &eye_x,
                               &eye_y,
                               &eye_z,
                               &eye_w);

  *x = eye_x;
  *y = eye_y;
  //*z = eye_z;
}

void
rut_camera_flush (RutCamera *camera)
{
  _rut_camera_flush_transforms (camera);

  if (camera->clear_fb)
    {
      cogl_framebuffer_clear4f (camera->fb,
                                COGL_BUFFER_BIT_COLOR |
                                COGL_BUFFER_BIT_DEPTH |
                                COGL_BUFFER_BIT_STENCIL,
                                cogl_color_get_red_float (&camera->bg_color),
                                cogl_color_get_green_float (&camera->bg_color),
                                cogl_color_get_blue_float (&camera->bg_color),
                                cogl_color_get_alpha_float (&camera->bg_color));
    }
}

void
rut_camera_end_frame (RutCamera *camera)
{
  double elapsed;
  CameraFlushState *state =
    cogl_object_get_user_data (COGL_OBJECT (camera->fb), &fb_camera_key);
  if (state)
    {
      if (state->in_frame != TRUE)
        g_warning ("Un-balanced rut_camera_flush/end frame calls. "
                   "_end before _flush");
      state->in_frame = FALSE;
    }
  camera->frame++;

  elapsed = g_timer_elapsed (camera->timer, NULL);
  if (elapsed > 1.0)
    {
      g_print ("fps = %f (camera = %p)\n",
               camera->frame / elapsed,
               camera);
      g_timer_start (camera->timer);
      camera->frame = 0;
    }
}

