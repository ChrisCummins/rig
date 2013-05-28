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

#include <config.h>

#include <math.h>

#include "rut-global.h"
#include "rut-types.h"
#include "rut-geometry.h"
#include "rut-mesh.h"
#include "rut-mesh-ply.h"

#include "components/rut-material.h"
#include "components/rut-model.h"

#define PI 3.14159265359

CoglPrimitive *
rut_model_get_primitive (RutObject *object)
{
  RutModel *model = object;

  if (!model->primitive)
    {
      if (model->mesh)
        {
          model->primitive =
            rut_mesh_create_primitive (model->ctx, model->mesh);
        }
    }

  return model->primitive;
}

RutType rut_model_type;

static RutComponentableVTable _rut_model_componentable_vtable = {
  .draw = NULL
};

static RutPrimableVTable _rut_model_primable_vtable = {
  .get_primitive = rut_model_get_primitive
};

static RutPickableVTable _rut_model_pickable_vtable = {
  .get_mesh = rut_model_get_mesh
};

static void
_rut_model_free (void *object)
{
  RutModel *model = object;

  if (model->primitive)
    cogl_object_unref (model->primitive);

  if (model->mesh)
    rut_refable_unref (model->mesh);

  g_slice_free (RutModel, model);
}

static RutRefCountableVTable _rut_model_ref_countable_vtable = {
  rut_refable_simple_ref,
  rut_refable_simple_unref,
  _rut_model_free
};

void
_rut_model_init_type (void)
{
  rut_type_init (&rut_model_type, "RigModel");
  rut_type_add_interface (&rut_model_type,
                          RUT_INTERFACE_ID_REF_COUNTABLE,
                          offsetof (RutModel, ref_count),
                          &_rut_model_ref_countable_vtable);
  rut_type_add_interface (&rut_model_type,
                          RUT_INTERFACE_ID_COMPONENTABLE,
                          offsetof (RutModel, component),
                          &_rut_model_componentable_vtable);
  rut_type_add_interface (&rut_model_type,
                          RUT_INTERFACE_ID_PRIMABLE,
                          0, /* no associated properties */
                          &_rut_model_primable_vtable);
  rut_type_add_interface (&rut_model_type,
                          RUT_INTERFACE_ID_PICKABLE,
                          0, /* no associated properties */
                          &_rut_model_pickable_vtable);
}

static RutModel *
_rut_model_new (RutContext *ctx)
{
  RutModel *model;

  model = g_slice_new0 (RutModel);
  rut_object_init (&model->_parent, &rut_model_type);
  model->ref_count = 1;
  model->component.type = RUT_COMPONENT_TYPE_GEOMETRY;
  model->ctx = ctx;

  return model;
}

static float
calculate_magnitude (float x,
                     float y,
                     float z)
{
  return sqrt (x * x + y * y + z * z);
}

float
calculate_dot_product (float *v1,
                       float *v2)
{
  return (v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2]);
}

float*
calculate_cross_product (float *v1,
                         float *v2)
{
  float *cross = g_new (float, 3);

  cross[0] = (v1[1] * v2[2]) - (v2[1] * v1[2]);
  cross[1] = (v1[2] * v2[0]) - (v2[2] * v1[0]);
  cross[2] = (v1[0] * v2[1]) - (v2[0] * v1[1]);

  return cross;
}

static void
normalize_vertex (float *vertex)
{
  float magnitude = calculate_magnitude (vertex[0], vertex[1], vertex[2]);

  vertex[0] = vertex[0] / magnitude;
  vertex[1] = vertex[1] / magnitude;
  vertex[2] = vertex[2] / magnitude;
}

static CoglBool
measure_mesh_x_cb (void **attribute_data,
                   int vertex_index,
                   void *user_data)
{
  RutModel *model = user_data;
  float *pos = attribute_data[0];

  if (pos[0] < model->min_x)
    model->min_x = pos[0];
  if (pos[0] > model->max_x)
    model->max_x = pos[0];

  return TRUE;
}

static CoglBool
measure_mesh_xy_cb (void **attribute_data,
                    int vertex_index,
                    void *user_data)
{
  RutModel *model = user_data;
  float *pos = attribute_data[0];

  measure_mesh_x_cb (attribute_data, vertex_index, user_data);

  if (pos[1] < model->min_y)
    model->min_y = pos[1];
  if (pos[1] > model->max_y)
    model->max_y = pos[1];

  return TRUE;
}

static CoglBool
measure_mesh_xyz_cb (void **attribute_data,
                     int vertex_index,
                     void *user_data)
{
  RutModel *model = user_data;
  float *pos = attribute_data[0];
  float *normal = attribute_data[1];
  float *tangent = attribute_data[2];

  measure_mesh_xy_cb (attribute_data, vertex_index, user_data);

  if (pos[2] < model->min_z)
    model->min_z = pos[2];
  if (pos[2] > model->max_z)
    model->max_z = pos[2];

  normal[0] = 0;
  normal[1] = 0;
  normal[2] = 0;

  tangent[0] = 0;
  tangent[1] = 0;
  tangent[2] = 0;

  return TRUE;
}

typedef struct _Vector
{
  float x, y, z;
  float nx, ny, nz;
  float tx, ty, tz;
  float s0, t0;
  float s1, t1;
}Vector;

typedef struct _InterVector
{
  Vector *original_vector;
  float flat_position[3];
  CoglBool uncovered;
}InterVector;

typedef struct _Polygon
{
  InterVector *vectors;
  int id;
  float tangent[3];
  float normal[3];
  float binormal[3];
}Polygon;

static void
generate_cylindrical_mapping (float *positions[3],
                              float *tex0[3],
                              float *tex1[3],
                              float *tex2[3],
                              float *tex3[3],
                              RutModel *model)
{
  float *vert_p0 = positions[0];
  float *vert_p1 = positions[1];
  float *vert_p2 = positions[2];
  float center[3], dir1[3], dir2[3], angle;
  int i;

  /* Vertex 1 */

  center[0] = (model->min_x + model->max_x) * 0.5;
  center[1] = vert_p0[1];
  center[2] = (model->min_z + model->max_z) * 0.5;

  dir2[0] = model->min_x - center[2];
  dir2[1] = vert_p0[1] - center[1];
  dir2[2] = model->min_z - center[2];

  dir1[0] = vert_p0[0] - center[0];
  dir1[1] = vert_p0[1] - center[1];
  dir1[2] = vert_p0[2] - center[2];

  angle = atan2 (dir1[0], dir1[2]) - atan2 (dir2[0], dir2[2]);

  if (angle < 0)
    angle = (2.0 * PI) + angle;

  if (angle > 0)
    tex0[0][0] = angle/ (2.0 * PI);
  else
    tex0[0][0] = 0;

  tex0[0][1] = (vert_p0[1] - model->min_y) / (model->max_y - model->min_y);

  if (tex0[0][0] <= 0)
    {
      if (vert_p0[0] < vert_p1[0] ||
        vert_p0[0] < vert_p2[0] )
        tex0[0][0] = 1;
    }

  if ((vert_p0[1] == model->max_y ||
       vert_p0[1] == model->min_y) &&
       vert_p0[1] == vert_p1[1] &&
       vert_p0[1] == vert_p2[1])
    {
      tex0[0][0] = (vert_p0[0] - model->min_x) /
      (model->max_x - model->min_x);
      tex0[0][1] = (vert_p0[2] - model->min_z) /
      (model->max_z - model->min_z);
    }

  /* Vertex 2 */

  center[1] = vert_p1[1];
  dir2[1] = vert_p1[1] - center[1];

  dir1[0] = vert_p1[0] - center[0];
  dir1[1] = vert_p1[1] - center[1];
  dir1[2] = vert_p1[2] - center[2];

  angle = atan2 (dir1[0], dir1[2]) - atan2 (dir2[0], dir2[2]);

  if (angle < 0)
    angle = (2.0 * PI) + angle;

  if (angle > 0)
    tex0[1][0] = angle/ (2.0 * PI);
  else
    tex0[1][0] = 0;

  tex0[1][1] = (vert_p1[1] - model->min_y) / (model->max_y - model->min_y);

  if (tex0[1][0] <= 0)
    {
      if (vert_p1[0] < vert_p0[0] ||
        vert_p1[0] < vert_p2[0] )
        tex0[1][0] = 1;
    }

    if ((vert_p1[1] == model->max_y ||
         vert_p1[1] == model->min_y) &&
         vert_p1[1] == vert_p0[1] &&
         vert_p1[1] == vert_p2[1])
    {
      tex0[1][0] = (vert_p1[0] - model->min_x) /
      (model->max_x - model->min_x);
      tex0[1][1] = (vert_p1[2] - model->min_z) /
      (model->max_z - model->min_z);
    }

  /* Vertex 3 */

  center[1] = vert_p2[1];
  dir2[1] = vert_p2[1] - center[1];

  dir1[0] = vert_p2[0] - center[0];
  dir1[1] = vert_p2[1] - center[1];
  dir1[2] = vert_p2[2] - center[2];

  angle = atan2 (dir1[0], dir1[2]) - atan2 (dir2[0], dir2[2]);

  if (angle < 0)
    angle = (2.0 * PI) + angle;

  if (angle > 0)
    tex0[2][0] = angle/ (2.0 * PI);
  else
    tex0[2][0] = 0;

  tex0[2][1] = (vert_p2[1] - model->min_y) / (model->max_y - model->min_y);

  if (tex0[2][0] <= 0)
    {
      if (vert_p2[0] < vert_p0[0] ||
        vert_p2[0] < vert_p1[0] )
        tex0[2][0] = 1;
    }

  if ((vert_p2[1] == model->max_y ||
       vert_p2[1] == model->min_y) &&
       vert_p2[1] == vert_p0[1] &&
       vert_p2[1] == vert_p1[1])
    {
      tex0[2][0] = (vert_p2[0] - model->min_x) /
      (model->max_x - model->min_x);
      tex0[2][1] = (vert_p2[2] - model->min_z) /
      (model->max_z - model->min_z);
    }

  for (i = 0; i < 3; i++)
    {
      tex1[i][0] = tex0[i][0];
      tex1[i][1] = tex0[i][1];
      tex2[i][0] = tex0[i][0];
      tex2[i][1] = tex0[i][1];
      tex3[i][0] = tex0[i][0];
      tex3[i][1] = tex0[i][1];
    }
}

static void
calculate_normals (float *positions[3],
                   float *normals[3],
                   Polygon *polygon)
{
  float *vert_p0 = positions[0];
  float *vert_p1 = positions[1];
  float *vert_p2 = positions[2];

  float *vert_n0 = normals[0];
  float *vert_n1 = normals[1];
  float *vert_n2 = normals[2];

  float edge1[3];
  float edge2[3];
  float *poly_normal;

  edge1[0] = vert_p1[0] - vert_p0[0];
  edge1[1] = vert_p1[1] - vert_p0[1];
  edge1[2] = vert_p1[2] - vert_p0[2];

  edge2[0] = vert_p2[0] - vert_p0[0];
  edge2[1] = vert_p2[1] - vert_p0[1];
  edge2[2] = vert_p2[2] - vert_p0[2];

  poly_normal = calculate_cross_product (edge1, edge2);

  normalize_vertex (poly_normal);

  vert_n0[0] += poly_normal[0];
  vert_n0[1] += poly_normal[1];
  vert_n0[2] += poly_normal[2];

  normalize_vertex (vert_n0);

  vert_n1[0] += poly_normal[0];
  vert_n1[1] += poly_normal[1];
  vert_n1[2] += poly_normal[2];

  normalize_vertex (vert_n1);

  vert_n2[0] += poly_normal[0];
  vert_n2[1] += poly_normal[1];
  vert_n2[2] += poly_normal[2];

  normalize_vertex (vert_n2);

  polygon->normal[0] = poly_normal[0];
  polygon->normal[1] = poly_normal[1];
  polygon->normal[2] = poly_normal[2];

  g_free (poly_normal);
}

static void
calculate_tangents (float *positions[3],
                    float *tangents[3],
                    float *tex[3],
                    Polygon *polygon)
{
  float *vert_p0 = positions[0];
  float *vert_p1 = positions[1];
  float *vert_p2 = positions[2];

  float *vert_t0 = tangents[0];
  float *vert_t1 = tangents[1];
  float *vert_t2 = tangents[2];

  float edge1[3];
  float edge2[3];
  float tex_edge1[2];
  float tex_edge2[2];
  float poly_tangent[3];
  float coef;

  edge1[0] = vert_p1[0] - vert_p0[0];
  edge1[1] = vert_p1[1] - vert_p0[1];
  edge1[2] = vert_p1[2] - vert_p0[2];

  edge2[0] = vert_p2[0] - vert_p0[0];
  edge2[1] = vert_p2[1] - vert_p0[1];
  edge2[2] = vert_p2[2] - vert_p0[2];

  tex_edge1[0] = tex[1][0] - tex[0][0];
  tex_edge1[1] = tex[1][1] - tex[0][1];

  tex_edge2[0] = tex[2][0] - tex[0][0];
  tex_edge2[1] = tex[2][1] - tex[0][1];

  coef = 1 / (tex_edge1[0] * tex_edge2[1] - tex_edge2[0] * tex_edge1[1]);

  poly_tangent[0] = coef * ((edge1[0] * tex_edge2[1]) -
  (edge2[0] * tex_edge1[1]));
  poly_tangent[1] = coef * ((edge1[1] * tex_edge2[1]) -
  (edge2[1] * tex_edge1[1]));
  poly_tangent[2] = coef * ((edge1[2] * tex_edge2[1]) -
  (edge2[2] * tex_edge1[1]));

  normalize_vertex (poly_tangent);

  vert_t0[0] += poly_tangent[0];
  vert_t0[1] += poly_tangent[1];
  vert_t0[2] += poly_tangent[2];

  normalize_vertex (vert_t0);

  vert_t1[0] += poly_tangent[0];
  vert_t1[1] += poly_tangent[1];
  vert_t1[2] += poly_tangent[2];

  normalize_vertex (vert_t1);

  vert_t2[0] += poly_tangent[0];
  vert_t2[1] += poly_tangent[1];
  vert_t2[2] += poly_tangent[2];

  normalize_vertex (vert_t2);

  polygon->tangent[0] = poly_tangent[0];
  polygon->tangent[1] = poly_tangent[1];
  polygon->tangent[2] = poly_tangent[2];
}

static CoglBool
calculate_tangent_space (void **attribute_data_v0,
                         void **attribute_data_v1,
                         void **attribute_data_v2,
                         int v0_index,
                         int v1_index,
                         int v2_index,
                         void *user_data)
{
  RutModel *model = user_data;
  Polygon *polygon = g_new (Polygon, 1);
  float *poly_binormal;

  float *positions[3] = { attribute_data_v0[0], attribute_data_v1[0],
                          attribute_data_v2[0] };
  float *normals[3] = { attribute_data_v0[1], attribute_data_v1[1],
                        attribute_data_v2[1] };
  float *tangents[3] = { attribute_data_v0[2], attribute_data_v1[2],
                         attribute_data_v2[2] };
  float *tex0[3] = { attribute_data_v0[3], attribute_data_v1[3],
                     attribute_data_v2[3] };
  float *tex1[3] = { attribute_data_v0[4], attribute_data_v1[4],
                     attribute_data_v2[4] };
  float *tex2[3] = { attribute_data_v0[5], attribute_data_v1[5],
                     attribute_data_v2[5] };
  float *tex3[3] = { attribute_data_v0[6], attribute_data_v1[6],
                     attribute_data_v2[6] };

  generate_cylindrical_mapping (positions, tex0, tex1, tex2, tex3, model);

  calculate_normals (positions, normals, polygon);

  calculate_tangents (positions, tangents, tex0, polygon);

  polygon->id = model->num_polygons;
  polygon->vectors = g_new (InterVector, 3);

  polygon->vectors[0].original_vector = g_new (Vector, 1);
  polygon->vectors[0].original_vector->x = positions[0][0];
  polygon->vectors[0].original_vector->y = positions[0][1];
  polygon->vectors[0].original_vector->z = positions[0][2];
  polygon->vectors[0].original_vector->nx = normals[0][0];
  polygon->vectors[0].original_vector->ny = normals[0][1];
  polygon->vectors[0].original_vector->nz = normals[0][2];
  polygon->vectors[0].original_vector->tx = tangents[0][0];
  polygon->vectors[0].original_vector->ty = tangents[0][1];
  polygon->vectors[0].original_vector->tz = tangents[0][2];
  polygon->vectors[0].original_vector->s0 = tex0[0][0];
  polygon->vectors[0].original_vector->t0 = tex0[0][1];
  polygon->vectors[0].uncovered = TRUE;

  polygon->vectors[1].original_vector = g_new (Vector, 1);
  polygon->vectors[1].original_vector->x = positions[1][0];
  polygon->vectors[1].original_vector->y = positions[1][1];
  polygon->vectors[1].original_vector->z = positions[1][2];
  polygon->vectors[1].original_vector->nx = normals[1][0];
  polygon->vectors[1].original_vector->ny = normals[1][1];
  polygon->vectors[1].original_vector->nz = normals[1][2];
  polygon->vectors[1].original_vector->tx = tangents[1][0];
  polygon->vectors[1].original_vector->ty = tangents[1][1];
  polygon->vectors[1].original_vector->tz = tangents[1][2];
  polygon->vectors[1].original_vector->s0 = tex0[1][0];
  polygon->vectors[1].original_vector->t0 = tex0[1][1];
  polygon->vectors[1].uncovered = TRUE;

  polygon->vectors[2].original_vector = g_new (Vector, 1);
  polygon->vectors[2].original_vector->x = positions[2][0];
  polygon->vectors[2].original_vector->y = positions[2][1];
  polygon->vectors[2].original_vector->z = positions[2][2];
  polygon->vectors[2].original_vector->nx = normals[2][0];
  polygon->vectors[2].original_vector->ny = normals[2][1];
  polygon->vectors[2].original_vector->nz = normals[2][2];
  polygon->vectors[2].original_vector->tx = tangents[2][0];
  polygon->vectors[2].original_vector->ty = tangents[2][1];
  polygon->vectors[2].original_vector->tz = tangents[2][2];
  polygon->vectors[2].original_vector->s0 = tex0[2][0];
  polygon->vectors[2].original_vector->t0 = tex0[2][1];
  polygon->vectors[2].uncovered = TRUE;

  poly_binormal = calculate_cross_product (polygon->normal,
                                           polygon->tangent);
  polygon->binormal[0] = poly_binormal[0];
  polygon->binormal[1] = poly_binormal[1];
  polygon->binormal[2] = poly_binormal[2];

  g_free (poly_binormal);

  g_queue_push_tail (model->polygons, (gpointer) polygon);
  model->num_polygons++;


  return TRUE;
}

RutModel *
rut_model_new_from_mesh (RutContext *ctx,
                         RutMesh *mesh)
{
  RutModel *model;
  RutAttribute *attribute;
  RutMeshVertexCallback measure_callback;

  model = _rut_model_new (ctx);
  model->type = RUT_MODEL_TYPE_FILE;
  model->mesh = rut_refable_ref (mesh);

  attribute = rut_mesh_find_attribute (model->mesh, "cogl_position_in");

  model->min_x = G_MAXFLOAT;
  model->max_x = G_MINFLOAT;
  model->min_y = G_MAXFLOAT;
  model->max_y = G_MINFLOAT;
  model->min_z = G_MAXFLOAT;
  model->max_z = G_MINFLOAT;

  if (attribute->n_components == 1)
    {
      model->min_y = model->max_y = 0;
      model->min_z = model->max_z = 0;
      measure_callback = measure_mesh_x_cb;
    }
  else if (attribute->n_components == 2)
    {
      model->min_z = model->max_z = 0;
      measure_callback = measure_mesh_xy_cb;
    }
  else if (attribute->n_components == 3)
    measure_callback = measure_mesh_xyz_cb;

  rut_mesh_foreach_vertex (model->mesh,
                           measure_callback,
                           model,
                           "cogl_position_in",
                           "cogl_normal_in",
                           "tangent_in",
                           NULL);

  return model;
}

typedef struct _AdjMatrixNode
{
  int is_connected;
  Polygon *origin;
  Polygon *destination;
}AdjMatrixNode;


static CoglBool
check_vertex_equality (Vector *v1, Vector *v2)
{
  return v1->x == v2->x &&
         v1->y == v2->y &&
         v1->z == v2->z;
}
static int
check_for_shared_edge (Polygon *poly1,
                       Polygon *poly2)
{
  int i, j;
  Vector *edges1[3][2];
  Vector *edges2[3][2];

  edges1[0][0] = poly1->vectors[0].original_vector;
  edges1[0][1] = poly1->vectors[1].original_vector;

  edges1[1][0] = poly1->vectors[1].original_vector;
  edges1[1][1] = poly1->vectors[2].original_vector;

  edges1[2][0] = poly1->vectors[2].original_vector;
  edges1[2][1] = poly1->vectors[0].original_vector;

  edges2[0][0] = poly2->vectors[0].original_vector;
  edges2[0][1] = poly2->vectors[1].original_vector;

  edges2[1][0] = poly2->vectors[1].original_vector;
  edges2[1][1] = poly2->vectors[2].original_vector;

  edges2[2][0] = poly2->vectors[2].original_vector;
  edges2[2][1] = poly2->vectors[0].original_vector;

  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
        {
          if (check_vertex_equality (edges1[i][0], edges2[j][0]) &&
              check_vertex_equality (edges1[i][1], edges2[j][1]))
            return 1;
            else if (check_vertex_equality (edges1[i][0], edges2[j][1]) &&
                     check_vertex_equality (edges1[i][1], edges2[j][0]))
            return 1;
        }
    }

  return 0;
}

static AdjMatrixNode**
generate_adjacency_matrix (RutModel *model)
{
  AdjMatrixNode **adj_matrix;
  int i, j;

  adj_matrix = g_new (AdjMatrixNode*, model->num_polygons);

  for (i = 0; i < model->num_polygons; i++)
    {
      Polygon *origin;
      adj_matrix[i] = g_new (AdjMatrixNode, model->num_polygons);
      origin = (Polygon*) g_queue_pop_head (model->polygons);
      g_queue_push_tail (model->polygons, origin);

      for (j = 0; j < model->num_polygons; j++)
        {
          Polygon *destination = (Polygon*) g_queue_pop_tail (model->polygons);

          adj_matrix[i][j].origin = origin;
          adj_matrix[i][j].destination = destination;
          if (origin->id == destination->id)
            adj_matrix[i][j].is_connected = 0;
          else
            adj_matrix[i][j].is_connected = check_for_shared_edge (origin,
                                                                   destination);
          g_queue_push_head (model->polygons, destination);
        }
    }

  return adj_matrix;
}

typedef struct _TexturePatch
{
  float width;
  float height;
  Polygon *center_polygon;
  Vector *center_vertex;
  GQueue *surface_patch;
}TexturePatch;

RutModel *
rut_model_new_from_asset (RutContext *ctx, RutAsset *asset)
{
  RutMesh *mesh = rut_asset_get_mesh (asset);
  RutModel *model;
  AdjMatrixNode **adj_matrix;
  int i;
  GQueue *nodes;
  Polygon *current_node;
  CoglBool *visited;
  TexturePatch *patch = g_new (TexturePatch, 1);

  if (!mesh)
    return NULL;

  model = rut_model_new_from_mesh (ctx, mesh);
  model->asset = rut_refable_ref (asset);

  model->polygons = g_queue_new ();
  nodes = g_queue_new ();
  model->num_polygons = 0;
  rut_mesh_foreach_triangle (model->mesh,
                             calculate_tangent_space,
                             model,
                             "cogl_position_in",
                             "cogl_normal_in",
                             "tangent_in",
                             "cogl_tex_coord0_in",
                             "cogl_tex_coord1_in",
                             "cogl_tex_coord2_in",
                             "cogl_tex_coord5_in",
                             NULL);

  adj_matrix = generate_adjacency_matrix (model);
  visited = g_new (CoglBool, model->num_polygons);

  for (i = 0; i < model->num_polygons; i++)
    visited[i] = FALSE;

  current_node = adj_matrix[0][0].origin;
  g_queue_push_tail (nodes, current_node);
  patch->width = 400;
  patch->height = 400;
  patch->center_polygon = current_node;
  patch->center_vertex = current_node->vectors[0].original_vector;
  patch->surface_patch = g_queue_new ();
  g_queue_push_tail (nodes, current_node);

  patch->center_polygon->vectors[0].flat_position[0] = patch->width / 2.0;
  patch->center_polygon->vectors[0].flat_position[1] = patch->height / 2.0;
  patch->center_polygon->vectors[0].original_vector->s1 = 0.5;
  patch->center_polygon->vectors[0].original_vector->t1 = 0.5;

  float edge1[3];
  float edge2[3];
  float s_displace[3];
  float t_displace[3];
  float magnitude;

  edge1[0] = patch->center_polygon->vectors[1].original_vector->x -
             patch->center_polygon->vectors[0].original_vector->x;

  edge1[1] = patch->center_polygon->vectors[1].original_vector->y -
             patch->center_polygon->vectors[0].original_vector->y;

  edge1[2] = patch->center_polygon->vectors[1].original_vector->z -
             patch->center_polygon->vectors[0].original_vector->z;

  s_displace[0] = edge1[0] * patch->center_polygon->binormal[0] +
                  

 /*
  magnitude = calculate_magnitude (edge1[0] *
                                   patch->center_polygon->binormal[0],
                                   edge1[1] *
                                   patch->center_polygon->binormal[1],
                                   edge1[2] *
                                   patch->center_polygon->binormal[2]);

  patch->center_polygon->vectors[1].flat_position[0] = magnitude +
    patch->center_polygon->vectors[0].flat_position[0];

  magnitude = calculate_magnitude (edge1[0] *
                                   patch->center_polygon->tangent[0],
                                   edge1[1] *
                                   patch->center_polygon->tangent[1],
                                   edge1[2] *
                                   patch->center_polygon->tangent[2]);

  patch->center_polygon->vectors[1].flat_position[1] = magnitude +
    patch->center_polygon->vectors[0].flat_position[1];

  patch->center_polygon->vectors[1].original_vector->s1 =
    patch->center_polygon->vectors[1].flat_position[0] / (float) patch->width;
  patch->center_polygon->vectors[1].original_vector->t1 =
    patch->center_polygon->vectors[1].flat_position[1] / (float) patch->height;*/

  edge2[0] = patch->center_polygon->vectors[2].original_vector->x -
             patch->center_polygon->vectors[0].original_vector->x;

  edge2[1] = patch->center_polygon->vectors[2].original_vector->y -
             patch->center_polygon->vectors[0].original_vector->y;

  edge2[2] = patch->center_polygon->vectors[2].original_vector->z -
             patch->center_polygon->vectors[0].original_vector->z;

  magnitude = calculate_magnitude (edge2[0] *
                                   patch->center_polygon->binormal[0],
                                   edge2[1] *
                                   patch->center_polygon->binormal[1],
                                   edge2[2] *
                                   patch->center_polygon->binormal[2]);

  patch->center_polygon->vectors[2].flat_position[0] = magnitude +
    patch->center_polygon->vectors[0].flat_position[0];

  magnitude = calculate_magnitude (edge2[0] *
                                   patch->center_polygon->tangent[0],
                                   edge2[1] *
                                   patch->center_polygon->tangent[1],
                                   edge2[2] *
                                   patch->center_polygon->tangent[2]);

  patch->center_polygon->vectors[2].flat_position[1] = magnitude +
    patch->center_polygon->vectors[0].flat_position[1];

  patch->center_polygon->vectors[2].original_vector->s1 =
    patch->center_polygon->vectors[2].flat_position[0] / (float) patch->width;
  patch->center_polygon->vectors[2].original_vector->t1 =
    patch->center_polygon->vectors[2].flat_position[1] / (float) patch->height;

  g_warning ("V1 %f %f\n", patch->center_polygon->vectors[0].flat_position[0],
    patch->center_polygon->vectors[0].flat_position[1]);
  g_warning ("V1 %f %f\n", patch->center_polygon->vectors[0].flat_position[0],
             patch->center_polygon->vectors[0].flat_position[1]);
  g_warning ("V2 %f %f\n", patch->center_polygon->vectors[1].flat_position[0],
             patch->center_polygon->vectors[1].flat_position[1]);
  g_warning ("V3 %f %f\n", patch->center_polygon->vectors[2].flat_position[0],
             patch->center_polygon->vectors[2].flat_position[1]);

/*
  while (!g_queue_is_empty (nodes))
    {
      current_node = g_queue_pop_head (nodes);
      if (!visited[current_node->id])
        {
          visited[current_node->id] = TRUE;
          g_queue_push_tail (patch->surface_patch, current_node);
          for (i = 0; i < model->num_polygons; i++)
            {
              if (adj_matrix[current_node->id][i].is_connected == 1)
                {
                  if (!visited[i])
                    {

                      g_queue_push_tail (nodes,
                        adj_matrix[current_node->id][i].destination);
                    }
                }
            }
        }
    }*/

  return model;
}

RutMesh *
rut_model_get_mesh (RutObject *self)
{
  RutModel *model = self;
  return model->mesh;
}

RutAsset *
rut_model_get_asset (RutModel *model)
{
  return model->asset;
}

