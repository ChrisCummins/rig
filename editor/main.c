#include <glib.h>
#include <gio/gio.h>
#include <sys/stat.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#include <cogl/cogl.h>

#include "rig.h"

//#define DEVICE_WIDTH 480.0
//#define DEVICE_HEIGHT 800.0
#define DEVICE_WIDTH 720.0
#define DEVICE_HEIGHT 1280.0

/*
 * Note: The size and padding for this circle texture have been carefully
 * chosen so it has a power of two size and we have enough padding to scale
 * down the circle to a size of 2 pixels and still have a 1 texel transparent
 * border which we rely on for anti-aliasing.
 */
#define CIRCLE_TEX_RADIUS 16
#define CIRCLE_TEX_PADDING 16

#define N_CUBES 5

#define INDENT_LEVEL 2

typedef struct _Data Data;

typedef struct _Node
{
  float t;
} Node;

typedef struct _NodeFloat
{
  float t;
  float value;
} NodeFloat;

typedef struct _NodeQuaternion
{
  float t;
  CoglQuaternion value;
} NodeQuaternion;

typedef struct _Path
{
  RigContext *ctx;
  RigProperty *progress_prop;
  RigProperty *prop;
  GQueue nodes;
  GList *pos;
} Path;

typedef struct _DiamondSlice
{
  RigObjectProps _parent;
  int ref_count;

  CoglMatrix rotate_matrix;

  CoglTexture *texture;

  float width;
  float height;

  CoglPipeline *pipeline;
  CoglPrimitive *primitive;

  RigGraphableProps graphable;
  RigPaintableProps paintable;

} DiamondSlice;

static uint8_t _diamond_slice_indices_data[] = {
    0,4,5,   0,5,1,  1,5,6,   1,6,2,   2,6,7,    2,7,3,
    4,8,9,   4,9,5,  5,9,10,  5,10,6,  6,10,11,  6,11,7,
    8,12,13, 8,13,9, 9,13,14, 9,14,10, 10,14,15, 10,15,11
};

typedef enum _AssetType {
  ASSET_TYPE_TEXTURE,
} AssetType;

#if 0
enum {
  ASSET_N_PROPS
};
#endif

typedef struct _Asset
{
  RigObjectProps _parent;

  Data *data;

  uint32_t id;

#if 0
  RigSimpleIntrospectableProps introspectable;
  RigProperty props[ASSET_N_PROPS];
#endif

  AssetType type;

  char *path;
  CoglTexture *texture;

} Asset;

#if 0
static RigPropertySpec _asset_prop_specs[] = {
  { 0 }
};
#endif

enum {
    ITEM_PROP_X,
    ITEM_PROP_Y,
    ITEM_PROP_Z,
    ITEM_PROP_ROTATION,
    ITEM_PROP_TRANSFORM,
    ITEM_PROP_RECEIVE_SHADOW,
    ITEM_PROP_CAST_SHADOW,
    ITEM_N_PROPS
};

typedef struct _Item
{
  RigObjectProps _parent;

  int ref_count;

  Data *data;

  uint32_t id;

  float x, y, z;
  CoglQuaternion rotation;
  CoglMatrix transform;

  Asset *texture_asset;

  CoglBool receive_shadow;
  CoglBool cast_shadow;

  CoglPipeline *pipeline;

  float diamond_size;
  DiamondSlice *diamond_slice;

  //float x0, y0, x1, y1;
  RigInputRegion *input_region;
  RigTransform *input_transform;

  RigGraphableProps graphable;
  RigPaintableProps paintable;

  RigProperty props[ITEM_N_PROPS];
  RigSimpleIntrospectableProps introspectable;

} Item;

static RigPropertySpec _item_prop_specs[] = {
  {
    .name = "x",
    .nick = "X",
    .type = RIG_PROPERTY_TYPE_FLOAT,
    .data_offset = offsetof (Item, x),
  },
  {
    .name = "y",
    .nick = "Y",
    .type = RIG_PROPERTY_TYPE_FLOAT,
    .data_offset = offsetof (Item, y)
  },
  {
    .name = "z",
    .nick = "Z",
    .type = RIG_PROPERTY_TYPE_FLOAT,
    .data_offset = offsetof (Item, z)
  },
  {
    .name = "rotation",
    .nick = "Rotation",
    .type = RIG_PROPERTY_TYPE_QUATERNION,
    .data_offset = offsetof (Item, rotation)
  },
  {
    .name = "transform",
    .nick = "Transform",
  },
  {
    .name = "receive-shadow",
    .nick = "Shadow",
    .type = RIG_PROPERTY_TYPE_BOOLEAN,
    .data_offset = offsetof (Item, receive_shadow)
  },
  {
    .name = "cast-shadow",
    .nick = "Cast Shadow",
    .type = RIG_PROPERTY_TYPE_BOOLEAN,
    .data_offset = offsetof (Item, cast_shadow)
  },
  { 0 }
};

enum {
  TRANSITION_PROP_PROGRESS,
  TRANSITION_N_PROPS
};

typedef struct _Transition
{
  RigObjectProps _parent;

  Data *data;

  uint32_t id;

  float progress;

  GList *paths;

  RigProperty props[TRANSITION_N_PROPS];
  RigSimpleIntrospectableProps introspectable;

} Transition;

static RigPropertySpec _transition_prop_specs[] = {
  {
    .name = "progress",
    .type = RIG_PROPERTY_TYPE_FLOAT,
    .data_offset = offsetof (Transition, progress)
  },
  { 0 }
};

typedef struct _TestPaintContext
{
  RigPaintContext _parent;

  Data *data;

  GList *camera_stack;

  gboolean shadow_pass;

} TestPaintContext;

#ifndef __ANDROID__

static gboolean _rig_handset_in_device_mode = FALSE;
static char **_rig_handset_remaining_args = NULL;

static const GOptionEntry rig_handset_entries[] =
{
  { "device-mode", 'd', 0, 0,
    &_rig_handset_in_device_mode, "Run in Device Mode" },
  { G_OPTION_REMAINING, 0, 0, G_OPTION_ARG_STRING_ARRAY,
    &_rig_handset_remaining_args, "Project" },
  { 0 }
};

static char *_rig_project_dir = NULL;

#endif /* __ANDROID__ */

static void
save (Data *data);

static void
load (Data *data, const char *file);

static Path *
transition_find_path (Transition *transition,
                      const char *property);

static void
path_lerp_property (Path *path, float t);

void
node_float_lerp (NodeFloat *a,
                 NodeFloat *b,
                 float t,
                 float *value)
{
  float range = b->t - a->t;
  float offset = t - a->t;
  float factor = offset / range;

  *value = a->value + (b->value - a->value) * factor;
}

void
node_quaternion_lerp (NodeQuaternion *a,
                      NodeQuaternion *b,
                      float t,
                      CoglQuaternion *value)
{
  float range = b->t - a->t;
  float offset = t - a->t;
  float factor = offset / range;

  cogl_quaternion_nlerp (value, &a->value, &b->value, factor);
}

void
node_free (void *node, void *user_data)
{
  RigPropertyType type = GPOINTER_TO_UINT (user_data);

  switch (type)
    {
      case RIG_PROPERTY_TYPE_FLOAT:
	g_slice_free (NodeFloat, node);
	break;

      case RIG_PROPERTY_TYPE_QUATERNION:
	g_slice_free (NodeQuaternion, node);
	break;

      default:
        g_warn_if_reached ();
    }
}

NodeFloat *
node_new_for_float (float t, float value)
{
  NodeFloat *node = g_slice_new (NodeFloat);
  node->t = t;
  node->value = value;
  return node;
}

NodeQuaternion *
node_new_for_quaternion (float t, float angle, float x, float y, float z)
{
  NodeQuaternion *node = g_slice_new (NodeQuaternion);
  node->t = t;
  cogl_quaternion_init (&node->value, angle, x, y, z);

  return node;
}

void
path_free (Path *path)
{
  g_queue_foreach (&path->nodes, node_free, GUINT_TO_POINTER (path->prop->spec->type));
  g_queue_clear (&path->nodes);
  rig_ref_countable_unref (path->ctx);
  g_slice_free (Path, path);
}


static void
update_path_property_cb (RigProperty *path_property, void *user_data)
{
  Path *path = user_data;
  float progress = rig_property_get_float (path->progress_prop);

  path_lerp_property (path, progress);
}

Path *
path_new_for_property (RigContext *ctx,
                       RigProperty *progress_prop,
                       RigProperty *path_prop)
{
  Path *path = g_slice_new (Path);

  path->ctx = rig_ref_countable_ref (ctx);

  path->progress_prop = progress_prop;
  path->prop = path_prop;

  g_queue_init (&path->nodes);
  path->pos = NULL;

  rig_property_set_binding (path_prop,
                            update_path_property_cb,
                            path,
                            progress_prop,
                            NULL);

  return path;
}

static GList *
nodes_find_less_than (GList *start, float t)
{
  GList *l;

  for (l = start; l; l = l->prev)
    {
      Node *node = l->data;
      if (node->t < t)
        return l;
    }

  return NULL;
}

static GList *
nodes_find_less_than_equal (GList *start, float t)
{
  GList *l;

  for (l = start; l; l = l->prev)
    {
      Node *node = l->data;
      if (node->t <= t)
        return l;
    }

  return NULL;
}

static GList *
nodes_find_greater_than (GList *start, float t)
{
  GList *l;

  for (l = start; l; l = l->next)
    {
      Node *node = l->data;
      if (node->t > t)
        return l;
    }

  return NULL;
}

static GList *
nodes_find_greater_than_equal (GList *start, float t)
{
  GList *l;

  for (l = start; l; l = l->next)
    {
      Node *node = l->data;
      if (node->t >= t)
        return l;
    }

  return NULL;
}

/* Finds 1 point either side of the given t using the direction to resolve
 * which points to choose if t corresponds to a specific node.
 */
static void
path_find_control_links2 (Path *path,
                          float t,
                          int direction,
                          GList **n0,
                          GList **n1)
{
  GList *pos;
  Node *pos_node;

  if (G_UNLIKELY (path->nodes.head == NULL))
    {
      *n0 = NULL;
      *n1= NULL;
      return;
    }

  if (G_UNLIKELY (path->pos == NULL))
    path->pos = path->nodes.head;

  pos = path->pos;
  pos_node = pos->data;

  /*
   * Note:
   *
   * A node with t exactly == t may only be considered as the first control
   * point moving in the current direction.
   */

  if (direction > 0)
    {
      if (pos_node->t > t)
        /* > --- T -------- PosT ---- */
        pos = nodes_find_less_than_equal (pos, t);
      else
        {
          /* > --- PosT -------- T ---- */
          GList *tmp = nodes_find_greater_than (pos, t);
          pos = tmp->prev;
        }

      *n0 = pos;
      *n1 = pos->next;
    }
  else
    {
      if (pos_node->t > t)
        {
          /* < --- T -------- PosT ---- */
          GList *tmp = nodes_find_less_than (pos, t);
          pos = tmp->next;
        }
      else
        /* < --- PosT -------- T ---- */
        pos = nodes_find_greater_than_equal (pos, t);

      *n0 = pos;
      *n1 = pos->prev;
    }

  path->pos = pos;
}

void
path_find_control_points2 (Path *path,
                           float t,
                           int direction,
                           GList **n0,
                           GList **n1)
{
  GList *l0, *l1;
  path_find_control_links2 (path, t, direction, &l0, &l1);
  *n0 = l0->data;
  *n1 = l1->data;
}

/* Finds 2 points either side of the given t using the direction to resolve
 * which points to choose if t corresponds to a specific node. */
void
path_find_control_points4 (Path *path,
                           float t,
                           int direction,
                           Node **n0,
                           Node **n1,
                           Node **n2,
                           Node **n3)
{
  GList *l1, *l2;

  path_find_control_links2 (path, t, direction, &l1, &l2);

  if (direction > 0)
    {
      *n0 = l1->prev->data;
      *n3 = l2->next->data;
    }
  else
    {
      *n0 = l1->next->data;
      *n3 = l2->prev->data;
    }

  *n1 = l1->data;
  *n2 = l2->data;
}

void
node_print (void *node, void *user_data)
{
  RigPropertyType type = GPOINTER_TO_UINT (user_data);
  switch (type)
    {
      case RIG_PROPERTY_TYPE_FLOAT:
	{
	  NodeFloat *f_node = (NodeFloat *)node;

	  g_print (" t = %f value = %f\n", f_node->t, f_node->value);
          break;
	}

      case RIG_PROPERTY_TYPE_QUATERNION:
	{
	  NodeQuaternion *q_node = (NodeQuaternion *)node;
	  const CoglQuaternion *q = &q_node->value;
	  g_print (" t = %f [%f (%f, %f, %f)]\n",
                   q_node->t,
                   q->w, q->x, q->y, q->z);
	  break;
	}

      default:
        g_warn_if_reached ();
    }
}

void
path_print (Path *path)
{
  g_print ("path=%p\n", path);
  g_queue_foreach (&path->nodes, node_print, GUINT_TO_POINTER (path->prop->spec->type));
}

static int
path_find_t_cb (gconstpointer a, gconstpointer b)
{
  const Node *node = a;
  const float *t = b;

  if (node->t == *t)
    return 0;

  return 1;
}

static int
path_node_sort_t_func (const Node *a,
                       const Node *b,
                       void *user_data)
{
  if (a->t == b->t)
    return 0;
  else if (a->t < b->t)
    return -1;
  else
    return 1;
}

void
path_insert_float (Path *path,
                   float t,
                   float value)
{
  GList *link;
  NodeFloat *node;

#if 0
  g_print ("BEFORE:\n");
  path_print (path);
#endif

  link = g_queue_find_custom (&path->nodes, &t, path_find_t_cb);

  if (link)
    {
      node = link->data;
      node->value = value;
    }
  else
    {
      node = node_new_for_float (t, value);
      g_queue_insert_sorted (&path->nodes, node,
                             (GCompareDataFunc)path_node_sort_t_func,
                             NULL);
    }

#if 0
  g_print ("AFTER:\n");
  path_print (path);
#endif
}

void
path_insert_quaternion (Path *path,
                        float t,
                        float angle,
                        float x,
                        float y,
                        float z)
{
  GList *link;
  NodeQuaternion *node;

#if 0
  g_print ("BEFORE:\n");
  path_print (path);
#endif

  link = g_queue_find_custom (&path->nodes, &t, path_find_t_cb);

  if (link)
    {
      node = link->data;
      cogl_quaternion_init (&node->value, angle, x, y, z);
    }
  else
    {
      node = node_new_for_quaternion (t, angle, x, y, z);
      g_queue_insert_sorted (&path->nodes, node, path_node_sort_t_func, NULL);
    }

#if 0
  g_print ("AFTER:\n");
  path_print (path);
#endif
}

static void
path_lerp_property (Path *path, float t)
{
  Node *n0, *n1;

  path_find_control_points2 (path, t, 1,
                             &n0,
                             &n1);

  switch (path->prop->spec->type)
    {
    case RIG_PROPERTY_TYPE_FLOAT:
      {
        float value;

        node_float_lerp (n0, n1, t, &value);
        rig_property_set_float (&path->ctx->property_ctx, path->prop, value);
        break;
      }
    case RIG_PROPERTY_TYPE_QUATERNION:
      {
        CoglQuaternion value;
        node_quaternion_lerp (n0, n1, t, &value);
        rig_property_set_quaternion (&path->ctx->property_ctx, path->prop, &value);
        break;
      }
    }
}

#if 0
typedef struct _CatmullRomNode
{
  float t;
  float point[3];
} CatmullRomNode;

typedef struct _CatmullRom
{
  GQueue nodes;
  GList *last_pos;
} CatmullRom;

CatmullRom *
catmull_rom_new (void)
{
  CatmullRom *cmr = g_slice_new (CatmullRom);

  g_queue_init (&cmr->nodes);
  cmr->last_pos = NULL;

  return cmr;
}

static void
free_node (CatmullRomNode *node)
{
  g_slice_free (CatmullRomNode, node);
}

void
catmull_rom_free (CatmullRom *cmr)
{
  g_queue_foreach (&cmr->nodes, free_node, NULL);
  g_slice_free (CatmullRom, cmr);
}

static int
catmull_rom_node_sort_func (const CatmullRomNode *a,
                            const CatmullRomNode *b,
                            void *user_data)
{
  if (a->t == b->t)
    return 0;
  else if (a->t < b->t)
    return -1;
  else
    return 1;
}

void
catmull_rom_add_point (CatmullRom *cmr,
                       float t,
                       float x,
                       float y,
                       float z)
{
  CatmullRomNode *node = g_slice_new (CatmullRomNode);
  node->t = t;
  node->point[0] = x;
  node->point[1] = y;
  node->point[2] = z;

  g_queue_insert_sorted (&cmr->nodes, node, catmull_rom_node_sort_func, NULL);
}

void
_catmull_rom_interpolate_points (float t,
                                 CatmullRomNode *n0,
                                 CatmullRomNode *n1,
                                 CatmullRomNode *n2,
                                 CatmullRomNode *n3,
                                 float *x,
                                 float *y,
                                 float *z)
{
  /* ref: http://www.codeproject.com/Articles/30838/Overhauser-Catmull-Rom-Splines-for-Camera-Animatio */

  float tt = t * t;
  float tt2 = 2.0 * tt;
  float ttt = tt * t;
  float ttt3 = 3.0 * ttt;

  float b1 = 0.5 * (-ttt + tt2 - t);
  float b2 = 0.5 * (ttt3 - 5.0 * tt + 2);
  float b3 = 0.5 * (-ttt3 + tt2 + tt2 + t);
  float b4 = 0.5 * (ttt - tt);

  *x = b1 * n0->point[0] + b2 * n1->point[0] + b3 * n2->point[0] + b4 * n3->point[0];
  *y = b1 * n0->point[1] + b2 * n1->point[1] + b3 * n2->point[1] + b4 * n3->point[1];
  *z = b1 * n0->point[2] + b2 * n1->point[2] + b3 * n2->point[2] + b4 * n3->point[2];
}

#if 0
void
catmull_rom_get_point (CatmullRom *cmr,
                       float t,
                       float *x,
                       float *y,
                       float *z)
{
  if (G_UNLIKELY (cmr->nodes->head == NULL))
    {
      *x = 0;
      *y = 0;
      *z = 0;
    }

  if (G_UNLIKELY (cmr->last_pos == NULL))
    cmr->last_pos = cmr->nodes->head;

  CatmullRomNode *pos = cmr->last_pos->head->data;

  if (pos->t == t)
  if (pos->t < t)
    {

    }

}
#endif

#endif

typedef struct _Bezier3D
{
  float p0[3];
  float p1[3];
  float p2[3];
  float p3[3];
} Bezier3D;

typedef enum _State
{
  STATE_NONE
} State;

enum {
  DATA_PROP_WIDTH,
  DATA_PROP_HEIGHT,
  //DATA_PROP_PATH_T,

  DATA_N_PROPS
};

struct _Data
{
  RigCamera *camera;
  RigObject *root;
  RigObject *scene;

  CoglPipeline *shadow_color_tex;
  CoglPipeline *shadow_map_tex;

  CoglPipeline *root_pipeline;
  CoglPipeline *default_pipeline;

  Bezier3D bezier;

  State state;

  RigShell *shell;
  RigContext *ctx;
  CoglOnscreen *onscreen;

  /* shadow mapping */
  CoglOffscreen *shadow_fb;
  CoglTexture2D *shadow_color;
  CoglTexture *shadow_map;
  RigCamera *shadow_map_camera;

  CoglIndices *diamond_slice_indices;
  CoglTexture *circle_texture;

  CoglTexture *light_icon;
  CoglTexture *clip_plane_icon;

  //float width;
  //RigProperty width_property;
  //float height;
  //RigProperty height_property;

  RigTransform *top_bar_transform;
  RigTransform *left_bar_transform;
  RigTransform *right_bar_transform;
  RigTransform *main_transform;
  RigTransform *bottom_bar_transform;

  RigTransform *screen_area_transform;

  CoglPrimitive *grid_prim;
  CoglAttribute *circle_node_attribute;
  int circle_node_n_verts;

  //RigTransform *slider_transform;
  //RigSlider *slider;
  //RigProperty *slider_progress;
  RigRectangle *rect;
  float width;
  float height;
  float top_bar_height;
  float left_bar_width;
  float right_bar_width;
  float bottom_bar_height;
  float grab_margin;
  float main_width;
  float main_height;
  float screen_area_width;
  float screen_area_height;

  RigRectangle *top_bar_rect;
  RigRectangle *left_bar_rect;
  RigRectangle *right_bar_rect;
  RigRectangle *bottom_bar_rect;

  RigUIViewport *assets_vp;
  RigGraph *assets_list;

  RigUIViewport *tool_vp;
  RigObject *tool_list;

  RigCamera *timeline_camera;
  RigInputRegion *timeline_input_region;
  float timeline_width;
  float timeline_height;
  float timeline_len;
  float timeline_scale;

  RigUIViewport *timeline_vp;

  float grab_timeline_vp_t;
  float grab_timeline_vp_y;

  CoglMatrix main_view;
  float z_2d;

  RigEntity *main_camera_rig0; /* move to origin */
  RigEntity *main_camera_rig1; /* armature rotate rotate */
  RigEntity *main_camera_rig2; /* negative offset */
  RigEntity *main_camera_rig3; /* armature length */

  RigEntity *main_camera;
  RigCamera *main_camera_component;
  float main_camera_z;
  RigInputRegion *main_input_region;

  RigEntity *plane;
  RigEntity *cubes[N_CUBES];
  RigEntity *light;

  RigArcball arcball;
  CoglQuaternion saved_rotation;
  float origin[3];
  float saved_origin[3];

  //RigTransform *screen_area_transform;
  RigTransform *device_transform;

  RigTimeline *timeline;
  RigProperty *timeline_elapsed;
  RigProperty *timeline_progress;

  float grab_x;
  float grab_y;
  RigInputCallback key_focus_callback;

  GList *assets;

  uint32_t items_next_id;
  uint32_t entity_next_id;
  GList *items;
  GList *entities;
  GList *lights;
  GList *pickables;
  GList *transitions;

  Item *selected_item;
  RigEntity *selected_entity;
  Transition *selected_transition;

  RigTool *tool;

  /* picking ray */
  CoglPipeline *picking_ray_color;
  CoglPrimitive *picking_ray;
  CoglBool debug_pick_ray;

  //Path *path;
  //float path_t;
  //RigProperty path_property;

  RigProperty properties[DATA_N_PROPS];

};

static RigPropertySpec data_propert_specs[] = {
  {
    .name = "width",
    .type = RIG_PROPERTY_TYPE_FLOAT,
    .data_offset = offsetof (Data, width)
  },
  {
    .name = "height",
    .type = RIG_PROPERTY_TYPE_FLOAT,
    .data_offset = offsetof (Data, height)
  },
#if 0
  {
    .name = "t",
    .type = RIG_PROPERTY_TYPE_FLOAT,
    .data_offset = offsetof (Data, path_t)
  },
#endif
  { 0 }
};

#if 1
static void
bezier_3d_get_point (Bezier3D *bezier,
                     float t,
                     float *result)
{
  float u = 1 - t;
  float tt = t * t;
  float uu = u * u;
  float uuu = uu * u;
  float ttt = tt * t;

  memcpy (result, bezier->p0, sizeof (float) * 3);
  cogl_vector3_multiply_scalar (result, uuu); /* first term */

#define ADD_TERM(R, FACTOR, P) \
  do { \
      float factor = (FACTOR); \
      R[0] += factor * P[0]; \
      R[1] += factor * P[1]; \
      R[2] += factor * P[2]; \
  } while (0)

  ADD_TERM (result, 3 * uu * t, bezier->p1); /* second term */
  ADD_TERM (result, 3 * u * tt, bezier->p2); /* third term */
  ADD_TERM (result, ttt, bezier->p3); /* fourth term */

#undef ADD_TERM
}
#endif

#if 0
XYZ CubicBezier(XYZ p0,XYZ p1,XYZ p2,XYZ p3, double mu)
{
   XYZ a,b,c,p;

   c.x = 3 * (p1.x - p0.x);
   c.y = 3 * (p1.y - p0.y);
   c.z = 3 * (p1.z - p0.z);
   b.x = 3 * (p2.x - p1.x) - c.x;
   b.y = 3 * (p2.y - p1.y) - c.y;
   b.z = 3 * (p2.z - p1.z) - c.z;
   a.x = p3.x - p0.x - c.x - b.x;
   a.y = p3.y - p0.y - c.y - b.y;
   a.z = p3.z - p0.z - c.z - b.z;

   p.x = a.x * mu * mu * mu + b.x * mu * mu + c.x * mu + p0.x;
   p.y = a.y * mu * mu * mu + b.y * mu * mu + c.y * mu + p0.y;
   p.z = a.z * mu * mu * mu + b.z * mu * mu + c.z * mu + p0.z;

   return(p);
}
#endif



static void
bezier_3d_init (Bezier3D *bezier,
                float *points)
{
  memcpy (bezier->p0, points, sizeof (float) * 3);
  memcpy (bezier->p1, points + 3, sizeof (float) * 3);
  memcpy (bezier->p2, points + 6, sizeof (float) * 3);
  memcpy (bezier->p3, points + 9, sizeof (float) * 3);
}

#if 0
#define BEZIER_MAX_RECURSE_DEPTH 16
static void
bezier_3d_sub (CoglPath *path,
               Bezier3D *cubic)
{
  Bezier3D cubics[BEZIER_MAX_RECURSE_DEPTH];
  Bezier3D *cleft;
  Bezier3D *cright;
  Bezier3D *c;
  float dif1[3];
  float dif2[3];
  float mm[3];
  float c1[3];
  float c2[3];
  float c3[3];
  float c4[3];
  float c5[3];
  int cindex;

  /* Put first curve on stack */
  cubics[0] = *cubic;
  cindex =  0;

  while (cindex >= 0)
    {
      c = &cubics[cindex];

      /* Calculate distance of control points from their
       * counterparts on the line between end points */
      dif1.x = (c->p2.x * 3) - (c->p1.x * 2) - c->p4.x;
      dif1.y = (c->p2.y * 3) - (c->p1.y * 2) - c->p4.y;
      dif2.x = (c->p3.x * 3) - (c->p4.x * 2) - c->p1.x;
      dif2.y = (c->p3.y * 3) - (c->p4.y * 2) - c->p1.y;

      if (dif1.x < 0)
        dif1.x = -dif1.x;
      if (dif1.y < 0)
        dif1.y = -dif1.y;
      if (dif2.x < 0)
        dif2.x = -dif2.x;
      if (dif2.y < 0)
        dif2.y = -dif2.y;


      /* Pick the greatest of two distances */
      if (dif1.x < dif2.x) dif1.x = dif2.x;
      if (dif1.y < dif2.y) dif1.y = dif2.y;

      /* Cancel if the curve is flat enough */
      if (dif1.x + dif1.y <= 1.0 ||
	  cindex == BEZIER_MAX_RECURSE_DEPTH-1)
	{
	  /* Add subdivision point (skip last) */
	  if (cindex == 0)
            return;

	  add_node (nodes, FALSE, c->p4.x, c->p4.y);

	  --cindex;

          continue;
	}

      /* Left recursion goes on top of stack! */
      cright = c; cleft = &cubics[++cindex];

      /* Subdivide into 2 sub-curves */
      c1.x = ((c->p1.x + c->p2.x) / 2);
      c1.y = ((c->p1.y + c->p2.y) / 2);
      mm.x = ((c->p2.x + c->p3.x) / 2);
      mm.y = ((c->p2.y + c->p3.y) / 2);
      c5.x = ((c->p3.x + c->p4.x) / 2);
      c5.y = ((c->p3.y + c->p4.y) / 2);

      c2.x = ((c1.x + mm.x) / 2);
      c2.y = ((c1.y + mm.y) / 2);
      c4.x = ((mm.x + c5.x) / 2);
      c4.y = ((mm.y + c5.y) / 2);

      c3.x = ((c2.x + c4.x) / 2);
      c3.y = ((c2.y + c4.y) / 2);

      /* Add left recursion to stack */
      cleft->p1 = c->p1;
      cleft->p2 = c1;
      cleft->p3 = c2;
      cleft->p4 = c3;

      /* Add right recursion to stack */
      cright->p1 = c3;
      cright->p2 = c4;
      cright->p3 = c5;
      cright->p4 = c->p4;
    }
}
#endif

#if 0
static UIViewport *
ui_viewport_new (float width,
                 float height,
                 float doc_x,
                 float doc_y,
                 float doc_x_scale,
                 float doc_y_scale)
{
  UIViewport *vp = g_slice_new (UIViewport);

  vp->width = width;
  vp->height = height;
  vp->doc_x = doc_x;
  vp->doc_y = doc_y;
  vp->doc_x_scale = doc_x_scale;
  vp->doc_y_scale = doc_y_scale;

  return vp;
}

static void
ui_viewport_free (UIViewport *vp)
{
  g_slice_free (UIViewport, vp);
}
#endif

static void
_diamond_slice_free (void *object)
{
  DiamondSlice *diamond_slice = object;

  cogl_object_unref (diamond_slice->texture);

  cogl_object_unref (diamond_slice->pipeline);
  cogl_object_unref (diamond_slice->primitive);

  g_slice_free (DiamondSlice, object);
}

RigRefCountableVTable _diamond_slice_ref_countable_vtable = {
  rig_ref_countable_simple_ref,
  rig_ref_countable_simple_unref,
  _diamond_slice_free
};

static RigGraphableVTable _diamond_slice_graphable_vtable = {
    0
};

static void
_diamond_slice_paint (RigObject *object,
                       RigPaintContext *paint_ctx)
{
  DiamondSlice *diamond_slice = object;
  RigCamera *camera = paint_ctx->camera;
  CoglFramebuffer *fb = rig_camera_get_framebuffer (camera);

  cogl_framebuffer_draw_primitive (fb,
                                   diamond_slice->pipeline,
                                   diamond_slice->primitive);
}

static RigPaintableVTable _diamond_slice_paintable_vtable = {
  _diamond_slice_paint
};

RigType _diamond_slice_type;

static void
_diamond_slice_init_type (void)
{
  rig_type_init (&_diamond_slice_type);
  rig_type_add_interface (&_diamond_slice_type,
                           RIG_INTERFACE_ID_REF_COUNTABLE,
                           offsetof (DiamondSlice, ref_count),
                           &_diamond_slice_ref_countable_vtable);
  rig_type_add_interface (&_diamond_slice_type,
                           RIG_INTERFACE_ID_GRAPHABLE,
                           offsetof (DiamondSlice, graphable),
                           &_diamond_slice_graphable_vtable);
  rig_type_add_interface (&_diamond_slice_type,
                           RIG_INTERFACE_ID_PAINTABLE,
                           offsetof (DiamondSlice, paintable),
                           &_diamond_slice_paintable_vtable);
}

typedef struct _VertexP2T2T2
{
  float x, y, s0, t0, s1, t1;
} VertexP2T2T2;

static CoglPrimitive *
primitive_new_p2t2t2 (CoglContext *ctx,
                      CoglVerticesMode mode,
                      int n_vertices,
                      const VertexP2T2T2 *data)
{
  CoglAttributeBuffer *attribute_buffer =
    cogl_attribute_buffer_new (ctx, n_vertices * sizeof (VertexP2T2T2), data);
  CoglAttribute *attributes[3];
  CoglPrimitive *primitive;
  int i;

  attributes[0] = cogl_attribute_new (attribute_buffer,
                                      "cogl_position_in",
                                      sizeof (VertexP2T2T2),
                                      offsetof (VertexP2T2T2, x),
                                      2,
                                      COGL_ATTRIBUTE_TYPE_FLOAT);
  attributes[1] = cogl_attribute_new (attribute_buffer,
                                      "cogl_tex_coord0_in",
                                      sizeof (VertexP2T2T2),
                                      offsetof (VertexP2T2T2, s0),
                                      2,
                                      COGL_ATTRIBUTE_TYPE_FLOAT);
  attributes[2] = cogl_attribute_new (attribute_buffer,
                                      "cogl_tex_coord1_in",
                                      sizeof (VertexP2T2T2),
                                      offsetof (VertexP2T2T2, s1),
                                      2,
                                      COGL_ATTRIBUTE_TYPE_FLOAT);

  cogl_object_unref (attribute_buffer);

  primitive = cogl_primitive_new_with_attributes (mode,
                                                  n_vertices,
                                                  attributes,
                                                  3);

  for (i = 0; i < 3; i++)
    cogl_object_unref (attributes[i]);

  return primitive;
}


DiamondSlice *
diamond_slice_new (Data *data,
                   CoglTexture *texture,
                   float size)
{
  RigContext *ctx = data->ctx;
  DiamondSlice *diamond_slice = g_slice_new (DiamondSlice);
  float width = size;
  float height = size;
#define DIAMOND_SLICE_CORNER_RADIUS 20
  CoglMatrix matrix;
  float tex_aspect;
  float t;
  float tex_width = cogl_texture_get_width (texture);
  float tex_height = cogl_texture_get_height (texture);

  rig_object_init (&diamond_slice->_parent, &_diamond_slice_type);

  diamond_slice->ref_count = 1;

  rig_graphable_init (RIG_OBJECT (diamond_slice));

  diamond_slice->texture = cogl_object_ref (texture);

  diamond_slice->width = width;
  diamond_slice->height = height;

  diamond_slice->pipeline = cogl_pipeline_new (ctx->cogl_context);
  cogl_pipeline_set_layer_texture (diamond_slice->pipeline, 0, data->circle_texture);
  cogl_pipeline_set_layer_texture (diamond_slice->pipeline, 1, texture);

    {
      //float tex_width = cogl_texture_get_width (rounded_texture);
      //float tex_height = cogl_texture_get_height (rounded_texture);

      /* x0,y0,x1,y1 and s0,t0,s1,t1 define the postion and texture
       * coordinates for the center rectangle... */
      float x0 = DIAMOND_SLICE_CORNER_RADIUS;
      float y0 = DIAMOND_SLICE_CORNER_RADIUS;
      float x1 = width - DIAMOND_SLICE_CORNER_RADIUS;
      float y1 = height - DIAMOND_SLICE_CORNER_RADIUS;

      /* The center region of the nine-slice can simply map to the
       * degenerate center of the circle */
      float s0 = 0.5;
      float t0 = 0.5;
      float s1 = 0.5;
      float t1 = 0.5;

#if 0
      float s0 = DIAMOND_SLICE_CORNER_RADIUS left / tex_width;
      float t0 = top / tex_height;
      float s1 = (tex_width - right) / tex_width;
      float t1 = (tex_height - bottom) / tex_height;
#endif

      int n_vertices;
      int i;

      /*
       * 0,0      x0,0      x1,0      width,0
       * 0,0      s0,0      s1,0      1,0
       * 0        1         2         3
       *
       * 0,y0     x0,y0     x1,y0     width,y0
       * 0,t0     s0,t0     s1,t0     1,t0
       * 4        5         6         7
       *
       * 0,y1     x0,y1     x1,y1     width,y1
       * 0,t1     s0,t1     s1,t1     1,t1
       * 8        9         10        11
       *
       * 0,height x0,height x1,height width,height
       * 0,1      s0,1      s1,1      1,1
       * 12       13        14        15
       */

      VertexP2T2T2 vertices[] =
        {
          { 0,  0, 0, 0, 0, 0 },
          { x0, 0, s0, 0, x0, 0},
          { x1, 0, s1, 0, x1, 0},
          { width, 0, 1, 0, width, 0},

          { 0, y0, 0, t0, 0, y0},
          { x0, y0, s0, t0, x0, y0},
          { x1, y0, s1, t0, x1, y0},
          { width, y0, 1, t0, width, y0},

          { 0, y1, 0, t1, 0, y1},
          { x0, y1, s0, t1, x0, y1},
          { x1, y1, s1, t1, x1, y1},
          { width, y1, 1, t1, width, y1},

          { 0, height, 0, 1, 0, height},
          { x0, height, s0, 1, x0, height},
          { x1, height, s1, 1, x1, height},
          { width, height, 1, 1, width, height},
        };

      cogl_matrix_init_identity (&diamond_slice->rotate_matrix);
      cogl_matrix_rotate (&diamond_slice->rotate_matrix, 45, 0, 0, 1);
      cogl_matrix_translate (&diamond_slice->rotate_matrix, - width / 2.0, - height / 2.0, 0);
      //cogl_matrix_translate (&diamond_slice->rotate_matrix, width / 2.0, height / 2.0, 0);

      n_vertices = sizeof (vertices) / sizeof (VertexP2T2T2);
      for (i = 0; i < n_vertices; i++)
        {
          float z = 0, w = 1;

          cogl_matrix_transform_point (&diamond_slice->rotate_matrix,
                                       &vertices[i].x,
                                       &vertices[i].y,
                                       &z,
                                       &w);
        }

      cogl_matrix_init_identity (&matrix);
      tex_aspect = tex_width / tex_height;

      t = 0.5 / sinf (G_PI_4);

      /* FIXME: hack */
      cogl_matrix_translate (&matrix, 0.5, 0, 0);

      if (tex_aspect < 1) /* taller than it is wide */
        {
          float s_scale = (1.0 / width) * t;

          float t_scale = s_scale * (1.0 / tex_aspect);

          cogl_matrix_scale (&matrix, s_scale, t_scale, t_scale);
        }
      else /* wider than it is tall */
        {
          float t_scale = (1.0 / height) * t;

          float s_scale = t_scale * tex_aspect;

          cogl_matrix_scale (&matrix, s_scale, t_scale, t_scale);
        }

      cogl_matrix_rotate (&matrix, 45, 0, 0, 1);

      n_vertices = sizeof (vertices) / sizeof (VertexP2T2T2);
      for (i = 0; i < n_vertices; i++)
        {
          float z = 0, w = 1;

          cogl_matrix_transform_point (&matrix,
                                       &vertices[i].s1,
                                       &vertices[i].t1,
                                       &z,
                                       &w);
        }


      diamond_slice->primitive =
        primitive_new_p2t2t2 (ctx->cogl_context,
                              COGL_VERTICES_MODE_TRIANGLES,
                              n_vertices,
                              vertices);

      /* The vertices uploaded only map to the key intersection points of the
       * 9-slice grid which isn't a topology that GPUs can handle directly so
       * this specifies an array of indices that allow the GPU to interpret the
       * vertices as a list of triangles... */
      cogl_primitive_set_indices (diamond_slice->primitive,
                                  data->diamond_slice_indices,
                                  sizeof (_diamond_slice_indices_data) /
                                  sizeof (_diamond_slice_indices_data[0]));
    }

  return diamond_slice;
}

CoglPrimitive *
create_grid (RigContext *ctx,
             float width,
             float height,
             float x_space,
             float y_space)
{
  GArray *lines = g_array_new (FALSE, FALSE, sizeof (CoglVertexP2));
  float x, y;
  int n_lines = 0;

  for (x = 0; x < width; x += x_space)
    {
      CoglVertexP2 p[2] = {
        { .x = x, .y = 0 },
        { .x = x, .y = height }
      };
      g_array_append_vals (lines, p, 2);
      n_lines++;
    }

  for (y = 0; y < height; y += y_space)
    {
      CoglVertexP2 p[2] = {
        { .x = 0, .y = y },
        { .x = width, .y = y }
      };
      g_array_append_vals (lines, p, 2);
      n_lines++;
    }

  return cogl_primitive_new_p2 (ctx->cogl_context,
                                COGL_VERTICES_MODE_LINES,
                                n_lines * 2,
                                (CoglVertexP2 *)lines->data);
}

static const float jitter_offsets[32] =
{
  0.375f, 0.4375f,
  0.625f, 0.0625f,
  0.875f, 0.1875f,
  0.125f, 0.0625f,

  0.375f, 0.6875f,
  0.875f, 0.4375f,
  0.625f, 0.5625f,
  0.375f, 0.9375f,

  0.625f, 0.3125f,
  0.125f, 0.5625f,
  0.125f, 0.8125f,
  0.375f, 0.1875f,

  0.875f, 0.9375f,
  0.875f, 0.6875f,
  0.125f, 0.3125f,
  0.625f, 0.8125f
};

/* XXX: This assumes that the primitive is being drawn in pixel coordinates,
 * since we jitter the modelview not the projection.
 */
static void
draw_jittered_primitive4f (Data *data,
                           CoglFramebuffer *fb,
                           CoglPrimitive *prim,
                           float red,
                           float green,
                           float blue)
{
  CoglPipeline *pipeline = cogl_pipeline_new (data->ctx->cogl_context);
  int i;

  cogl_pipeline_set_color4f (pipeline,
                             red / 16.0f,
                             green / 16.0f,
                             blue / 16.0f,
                             1.0f / 16.0f);

  for (i = 0; i < 16; i++)
    {
      const float *offset = jitter_offsets + 2 * i;

      cogl_framebuffer_push_matrix (fb);
      cogl_framebuffer_translate (fb, offset[0], offset[1], 0);
      cogl_framebuffer_draw_primitive (fb, pipeline, prim);
      cogl_framebuffer_pop_matrix (fb);
    }

  cogl_object_unref (pipeline);
}

static void
camera_update_view (Data *data, RigEntity *camera, CoglBool shadow_map)
{
  RigCamera *camera_component =
    rig_entity_get_component (camera, RIG_COMPONENT_TYPE_CAMERA);
  CoglMatrix transform;
  CoglMatrix inverse_transform;
  CoglMatrix view;
  CoglMatrix tmp;

  /* translate to z_2d and scale */
  view = data->main_view;

  /* translate into main_area */
  cogl_matrix_multiply (&view, &view, rig_transform_get_matrix (data->screen_area_transform));

  /* scale to device coordinates */
  cogl_matrix_multiply (&view, &view, rig_transform_get_matrix (data->device_transform));

#if 1
  rig_graphable_get_transform (camera, &transform);
  cogl_matrix_get_inverse (&transform, &inverse_transform);

  //tmp = view;
  //cogl_matrix_multiply (&view, &inverse_transform, &tmp);

  /* apply the camera viewing transform */
  cogl_matrix_multiply (&view, &view, &inverse_transform);
#endif

  //view = inverse_transform;

  if (shadow_map)
    {
      CoglMatrix flipped_view;
      cogl_matrix_init_identity (&flipped_view);
      cogl_matrix_scale (&flipped_view, 1, -1, 1);
      cogl_matrix_multiply (&flipped_view, &flipped_view, &view);
      rig_camera_set_view_transform (camera_component, &flipped_view);
    }
  else
    {
      CoglMatrix tmp = view;
      //cogl_matrix_multiply (&view, &view, &data->main_view);
      //cogl_matrix_multiply (&view, &data->main_view, &tmp);
      rig_camera_set_view_transform (camera_component, &view);
    }
}

static void
get_normal_matrix (const CoglMatrix *matrix,
                   float *normal_matrix)
{
  CoglMatrix inverse_matrix;

  /* Invert the matrix */
  cogl_matrix_get_inverse (matrix, &inverse_matrix);

  /* Transpose it while converting it to 3x3 */
  normal_matrix[0] = inverse_matrix.xx;
  normal_matrix[1] = inverse_matrix.xy;
  normal_matrix[2] = inverse_matrix.xz;

  normal_matrix[3] = inverse_matrix.yx;
  normal_matrix[4] = inverse_matrix.yy;
  normal_matrix[5] = inverse_matrix.yz;

  normal_matrix[6] = inverse_matrix.zx;
  normal_matrix[7] = inverse_matrix.zy;
  normal_matrix[8] = inverse_matrix.zz;
}

CoglPipeline *
get_entity_pipeline (Data *data,
                     RigEntity *entity,
                     RigComponent *geometry,
                     gboolean shadow_pass)
{
  CoglSnippet *snippet;
  CoglDepthState depth_state;
  RigMaterial *material =
    rig_entity_get_component (entity, RIG_COMPONENT_TYPE_MATERIAL);
  CoglPipeline *pipeline;
  GList *l;

  pipeline = rig_entity_get_pipeline_cache (entity);
  if (pipeline)
    return cogl_object_ref (pipeline);

  pipeline = cogl_pipeline_new (data->ctx->cogl_context);

#if 0
  /* NB: Our texture colours aren't premultiplied */
  cogl_pipeline_set_blend (pipeline,
                           "RGB = ADD(SRC_COLOR*(SRC_COLOR[A]), DST_COLOR*(1-SRC_COLOR[A]))"
                           "A   = ADD(SRC_COLOR, DST_COLOR*(1-SRC_COLOR[A]))",
                           NULL);
#endif

#if 0
  if (rig_object_get_type (geometry) == &rig_diamond_type)
    rig_geometry_component_update_pipeline (geometry, pipeline);

  for (l = data->lights; l; l = l->next)
    light_update_pipeline (l->data, pipeline);

  pipeline = cogl_pipeline_new (rig_cogl_context);
#endif

  cogl_pipeline_set_color4f (pipeline, 0.8f, 0.8f, 0.8f, 1.f);

  /* enable depth testing */
  cogl_depth_state_init (&depth_state);
  cogl_depth_state_set_test_enabled (&depth_state, TRUE);
  cogl_pipeline_set_depth_state (pipeline, &depth_state, NULL);

  /* Vertex shader setup for lighting */
  snippet = cogl_snippet_new (COGL_SNIPPET_HOOK_VERTEX,

      /* definitions */
      "uniform mat3 normal_matrix;\n"
      "varying vec3 normal_direction, eye_direction;\n",

      /* post */
      "normal_direction = normalize(normal_matrix * cogl_normal_in);\n"
      //"normal_direction = cogl_normal_in;\n"
      "eye_direction    = -vec3(cogl_modelview_matrix * cogl_position_in);\n"
  );

  cogl_pipeline_add_snippet (pipeline, snippet);
  cogl_object_unref (snippet);

#if 0
  /* Vertex shader setup for shadow mapping */
  snippet = cogl_snippet_new (COGL_SNIPPET_HOOK_VERTEX,

      /* definitions */
      "uniform mat4 light_shadow_matrix;\n"
      "varying vec4 shadow_coords;\n",

      /* post */
      "shadow_coords = light_shadow_matrix * cogl_modelview_matrix *\n"
      "                cogl_position_in;\n"
  );

  cogl_pipeline_add_snippet (pipeline, snippet);
  cogl_object_unref (snippet);
#endif

  /* and fragment shader */
  snippet = cogl_snippet_new (COGL_SNIPPET_HOOK_FRAGMENT,
      /* definitions */
      //"varying vec3 normal_direction;\n",
      "varying vec3 normal_direction, eye_direction;\n",
      /* post */
      "");
  //cogl_snippet_set_pre (snippet, "cogl_color_out = cogl_color_in;\n");

  cogl_pipeline_add_snippet (pipeline, snippet);
  cogl_object_unref (snippet);

  snippet = cogl_snippet_new (COGL_SNIPPET_HOOK_FRAGMENT,
      /* definitions */
      "uniform vec4 light0_ambient, light0_diffuse, light0_specular;\n"
      "uniform vec3 light0_direction_norm;\n",

      /* post */
      "vec4 final_color;\n"

      "vec3 L = light0_direction_norm;\n"
      "vec3 N = normalize(normal_direction);\n"

      "if (cogl_color_out.a <= 0.0)\n"
      "  discard;\n"

      "final_color = light0_ambient * cogl_color_out;\n"
      "float lambert = dot(N, L);\n"
      //"float lambert = 1.0;\n"

      "if (lambert > 0.0)\n"
      "{\n"
      "  final_color += cogl_color_out * light0_diffuse * lambert;\n"
      //"  final_color +=  vec4(1.0, 0.0, 0.0, 1.0) * light0_diffuse * lambert;\n"

      "  vec3 E = normalize(eye_direction);\n"
      "  vec3 R = reflect (-L, N);\n"
      "  float specular = pow (max(dot(R, E), 0.0),\n"
      "                        2.);\n"
      "  final_color += light0_specular * vec4(.6, .6, .6, 1.0) * specular;\n"
      "}\n"

      "cogl_color_out = final_color;\n"
      //"cogl_color_out = vec4(1.0, 0.0, 0.0, 1.0);\n"
  );

  cogl_pipeline_add_snippet (pipeline, snippet);
  cogl_object_unref (snippet);


#if 0

  /* Handle shadow mapping */

  snippet = cogl_snippet_new (COGL_SNIPPET_HOOK_FRAGMENT,
      /* declarations */
      "",

      /* post */
      "shadow_coords_d = shadow_coords / shadow_coords.w;\n"
      "cogl_texel7 =  cogl_texture_lookup7 (cogl_sampler7, cogl_tex_coord_in[0]);\n"
      "float distance_from_light = cogl_texel7.z + 0.0005;\n"
      "float shadow = 1.0;\n"
      "if (shadow_coords.w > 0.0 && distance_from_light < shadow_coords_d.z)\n"
      "    shadow = 0.5;\n"

      "cogl_color_out = shadow * cogl_color_out;\n"
  );

  cogl_pipeline_add_snippet (pipeline, snippet);
  cogl_object_unref (snippet);

  /* Hook the shadow map sampling */

  cogl_pipeline_set_layer_texture (pipeline, 7, data->shadow_map);

  snippet = cogl_snippet_new (COGL_SNIPPET_HOOK_TEXTURE_LOOKUP,
                              /* declarations */
                              "varying vec4 shadow_coords;\n"
                              "vec4 shadow_coords_d;\n",
                              /* post */
                              "");

  cogl_snippet_set_replace (snippet,
                            "cogl_texel = texture2D(cogl_sampler7, shadow_coords_d.st);\n");

  cogl_pipeline_add_layer_snippet (pipeline, 7, snippet);
  cogl_object_unref (snippet);
#endif

#if 1
  {
    RigLight *light = rig_entity_get_component (data->light, RIG_COMPONENT_TYPE_LIGHT);
    rig_light_set_uniforms (light, pipeline);
  }
#endif

#if 1
  if (rig_object_get_type (geometry) == &rig_diamond_type)
    {
      //pipeline = cogl_pipeline_new (data->ctx->cogl_context);

      //cogl_pipeline_set_depth_state (pipeline, &depth_state, NULL);

      rig_diamond_apply_mask (geometry, pipeline);

      //cogl_pipeline_set_color4f (pipeline, 1, 0, 0, 1);

#if 1
      if (material)
        {
          CoglPipeline *tmp = rig_material_get_pipeline (material);
          CoglTexture *tex = cogl_pipeline_get_layer_texture (tmp, 0);
          if (tex)
            {
              cogl_pipeline_set_layer_texture (pipeline, 1, tex);
            }
        }
#endif
    }
#endif

  rig_entity_set_pipeline_cache (entity, pipeline);

  return pipeline;
}

static RigTraverseVisitFlags
_rig_entitygraph_pre_paint_cb (RigObject *object,
                               int depth,
                               void *user_data)
{
  TestPaintContext *test_paint_ctx = user_data;
  RigPaintContext *paint_ctx = user_data;
  RigCamera *camera = paint_ctx->camera;
  CoglFramebuffer *fb = rig_camera_get_framebuffer (camera);
  RigPaintableVTable *vtable =
    rig_object_get_vtable (object, RIG_INTERFACE_ID_PAINTABLE);

  if (rig_object_is (object, RIG_INTERFACE_ID_TRANSFORMABLE))
    {
      const CoglMatrix *matrix = rig_transformable_get_matrix (object);
      cogl_framebuffer_push_matrix (fb);
      cogl_framebuffer_transform (fb, matrix);
    }

  if (rig_object_get_type (object) == &rig_entity_type)
    {
      RigComponent *geometry =
        rig_entity_get_component (object, RIG_COMPONENT_TYPE_GEOMETRY);
      CoglPipeline *pipeline;
      CoglPrimitive *primitive;
      CoglMatrix modelview_matrix;
      float normal_matrix[9];

      if (!geometry)
        {
          rig_entity_draw (object, fb);
          return RIG_TRAVERSE_VISIT_CONTINUE;
        }
#if 1
      pipeline = get_entity_pipeline (test_paint_ctx->data,
                                      object,
                                      geometry,
                                      test_paint_ctx->shadow_pass);
#endif

      primitive = rig_primable_get_primitive (geometry);

#if 1
      cogl_framebuffer_get_modelview_matrix (fb, &modelview_matrix);
      get_normal_matrix (&modelview_matrix, normal_matrix);

      {
        int location = cogl_pipeline_get_uniform_location (pipeline, "normal_matrix");
        cogl_pipeline_set_uniform_matrix (pipeline,
                                          location,
                                          3, /* dimensions */
                                          1, /* count */
                                          FALSE, /* don't transpose again */
                                          normal_matrix);
      }
#endif

      cogl_framebuffer_draw_primitive (fb,
                                       pipeline,
                                       primitive);

      /* FIXME: cache the pipeline with the entity */
      cogl_object_unref (pipeline);

#if 0
      geometry = rig_entity_get_component (object, RIG_COMPONENT_TYPE_GEOMETRY);
      material = rig_entity_get_component (object, RIG_COMPONENT_TYPE_MATERIAL);
      if (geometry && material)
        {
          if (rig_object_get_type (geometry) == &rig_diamond_type)
            {
              TestPaintContext *test_paint_ctx = paint_ctx;
              Data *data = test_paint_ctx->data;
              RigDiamondSlice *slice = rig_diamond_get_slice (geometry);
              CoglPipeline *template = rig_diamond_slice_get_pipeline_template (slice);
              CoglPipeline *material_pipeline = rig_material_get_pipeline (material);
              CoglPipeline *pipeline = cogl_pipeline_copy (template);
              //CoglPipeline *pipeline = cogl_pipeline_copy (data->root_pipeline);
              //CoglPipeline *pipeline = cogl_pipeline_new (data->ctx->cogl_context);

              /* FIXME: we should be combining the material and
               * diamond slice state together before now! */
              cogl_pipeline_set_layer_texture (pipeline, 1,
                                               cogl_pipeline_get_layer_texture (material_pipeline, 0));

              cogl_framebuffer_draw_primitive (fb,
                                               pipeline,
                                               slice->primitive);

              cogl_object_unref (pipeline);
            }
        }
#endif
      return RIG_TRAVERSE_VISIT_CONTINUE;
    }

  /* XXX:
   * How can we maintain state between the pre and post stages?  Is it
   * ok to just "sub-class" the paint context and maintain a stack of
   * state that needs to be shared with the post paint code.
   */

  return RIG_TRAVERSE_VISIT_CONTINUE;
}

static RigTraverseVisitFlags
_rig_entitygraph_post_paint_cb (RigObject *object,
                                int depth,
                                void *user_data)
{
  if (rig_object_is (object, RIG_INTERFACE_ID_TRANSFORMABLE))
    {
      RigPaintContext *paint_ctx = user_data;
      CoglFramebuffer *fb = rig_camera_get_framebuffer (paint_ctx->camera);
      cogl_framebuffer_pop_matrix (fb);
    }

  return RIG_TRAVERSE_VISIT_CONTINUE;
}

static void
compute_light_shadow_matrix (Data       *data,
                             CoglMatrix *light_matrix,
                             CoglMatrix *light_projection,
                             RigEntity  *light)
{
  CoglMatrix *main_camera, *light_transform, light_view;
  /* Move the unit data from [-1,1] to [0,1], column major order */
  float bias[16] = {
    .5f, .0f, .0f, .0f,
    .0f, .5f, .0f, .0f,
    .0f, .0f, .5f, .0f,
    .5f, .5f, .5f, 1.f
  };

  main_camera = rig_entity_get_transform (data->main_camera);
  light_transform = rig_entity_get_transform (light);
  cogl_matrix_get_inverse (light_transform, &light_view);

  cogl_matrix_init_from_array (light_matrix, bias);
  cogl_matrix_multiply (light_matrix, light_matrix, light_projection);
  cogl_matrix_multiply (light_matrix, light_matrix, &light_view);
  cogl_matrix_multiply (light_matrix, light_matrix, main_camera);
}

#if 1
static void
paint_main_area_camera (RigEntity *camera, TestPaintContext *test_paint_ctx)
{
  RigCamera *camera_component =
    rig_entity_get_component (camera, RIG_COMPONENT_TYPE_CAMERA);
  RigPaintContext *paint_ctx = &test_paint_ctx->_parent;
  Data *data = test_paint_ctx->data;
  CoglContext *ctx = data->ctx->cogl_context;
  CoglFramebuffer *fb = rig_camera_get_framebuffer (camera_component);
  GList *l;
  RigComponent *light;
  CoglFramebuffer *shadow_fb;

  camera_update_view (data, camera, FALSE);

  rig_camera_flush (camera_component);

  light = rig_entity_get_component (camera, RIG_COMPONENT_TYPE_LIGHT);
  test_paint_ctx->shadow_pass = light ? TRUE : FALSE;

#if 0
  {
    CoglPipeline *pipeline = cogl_pipeline_new (data->ctx->cogl_context);
    CoglMatrix view;
    float fovy = 10; /* y-axis field of view */
    float aspect = (float)data->main_width/(float)data->main_height;
    float z_near = 10; /* distance to near clipping plane */
    float z_far = 100; /* distance to far clipping plane */
#if 1
    fovy = 60;
    z_near = 1.1;
    z_far = 100;
#endif

    g_assert (data->main_width == cogl_framebuffer_get_viewport_width (fb));
    g_assert (data->main_height == cogl_framebuffer_get_viewport_height (fb));

    cogl_matrix_init_identity (&view);
    cogl_matrix_view_2d_in_perspective (&view,
                                        fovy, aspect, z_near, data->z_2d,
                                        DEVICE_WIDTH,
                                        DEVICE_HEIGHT);
    cogl_framebuffer_set_modelview_matrix (fb, &view);

    cogl_framebuffer_draw_rectangle (fb, pipeline,
                                     1, 1, DEVICE_WIDTH - 2, DEVICE_HEIGHT - 2);
    cogl_object_unref (pipeline);
  }
#endif

#if 0
  cogl_framebuffer_transform (fb, rig_transform_get_matrix (data->screen_area_transform));
  cogl_framebuffer_transform (fb, rig_transform_get_matrix (data->device_transform));
#endif

  if (!test_paint_ctx->shadow_pass)
    {
      CoglPipeline *pipeline = cogl_pipeline_new (ctx);
      cogl_pipeline_set_color4f (pipeline, 0, 0, 0, 1.0);
      cogl_framebuffer_draw_rectangle (fb,
                                       pipeline,
                                       0, 0, DEVICE_WIDTH, DEVICE_HEIGHT);
                                       //0, 0, data->pane_width, data->pane_height);
      cogl_object_unref (pipeline);

#if 1
      for (l = data->items; l; l = l->next)
        {
          Item *item = l->data;
          const CoglMatrix *matrix = &item->transform;

          cogl_framebuffer_push_matrix (fb);
          cogl_framebuffer_transform (fb, matrix);

          rig_paintable_paint (item->diamond_slice, paint_ctx);

          cogl_framebuffer_pop_matrix (fb);
        }
#endif
    }

  shadow_fb = COGL_FRAMEBUFFER (data->shadow_fb);

  /* update uniforms in materials */
  {
    CoglMatrix light_shadow_matrix, light_projection;
    CoglPipeline *pipeline;
    RigMaterial *material;
    const float *light_matrix;
    int location;

    cogl_framebuffer_get_projection_matrix (shadow_fb, &light_projection);
    compute_light_shadow_matrix (data,
                                 &light_shadow_matrix,
                                 &light_projection,
                                 data->light);
    light_matrix = cogl_matrix_get_array (&light_shadow_matrix);

    /* plane material */
    material = rig_entity_get_component (data->plane,
                                         RIG_COMPONENT_TYPE_MATERIAL);
    pipeline = rig_material_get_pipeline (material);
    location = cogl_pipeline_get_uniform_location (pipeline,
                                                   "light_shadow_matrix");
    cogl_pipeline_set_uniform_matrix (pipeline,
                                      location,
                                      4, 1,
                                      FALSE,
                                      light_matrix);

    /* cubes material */
    material = rig_entity_get_component (data->cubes[0],
                                         RIG_COMPONENT_TYPE_MATERIAL);
    pipeline = rig_material_get_pipeline (material);
    location = cogl_pipeline_get_uniform_location (pipeline,
                                                   "light_shadow_matrix");
    cogl_pipeline_set_uniform_matrix (pipeline,
                                      location,
                                      4, 1,
                                      FALSE,
                                      light_matrix);
  }



  rig_graphable_traverse (data->scene,
                          RIG_TRAVERSE_DEPTH_FIRST,
                          _rig_entitygraph_pre_paint_cb,
                          _rig_entitygraph_post_paint_cb,
                          test_paint_ctx);

  if (!test_paint_ctx->shadow_pass)
    {
      if (data->debug_pick_ray && data->picking_ray)
      //if (data->picking_ray)
        {
          cogl_framebuffer_draw_primitive (fb,
                                           data->picking_ray_color,
                                           data->picking_ray);
        }
#if 0
      for (l = data->entities; l; l = l->next)
        {
          RigEntity *entity = l->data;
          const CoglMatrix *transform;
          cogl_framebuffer_push_matrix (fb);

          transform = rig_entity_get_transform (entity);
          cogl_framebuffer_transform (fb, transform);

          rig_entity_draw (entity, fb);

          cogl_framebuffer_pop_matrix (fb);
        }
#endif

      draw_jittered_primitive4f (data, fb, data->grid_prim, 0.5, 0.5, 0.5);

      if (data->selected_entity)
        {
          rig_tool_update (data->tool, data->selected_entity);
          rig_tool_draw (data->tool, fb);
        }
    }

  rig_camera_end_frame (camera_component);
}
#endif

static void
paint_timeline_camera (RigCamera *camera, void *user_data)
{
  Data *data = user_data;
  CoglContext *ctx = data->ctx->cogl_context;
  CoglFramebuffer *fb = rig_camera_get_framebuffer (camera);
  GList *l;

  rig_camera_flush (camera);

  //cogl_framebuffer_push_matrix (fb);
  //cogl_framebuffer_transform (fb, rig_transformable_get_matrix (camera));

  if (data->selected_item)
    {
      //CoglContext *ctx = data->ctx->cogl_context;
      Item *item = data->selected_item;
      //int i;

      float viewport_x = 0;
      float viewport_y = 0;

      float viewport_t_scale =
        rig_ui_viewport_get_doc_scale_x (data->timeline_vp) *
        data->timeline_scale;

      float viewport_y_scale =
        rig_ui_viewport_get_doc_scale_y (data->timeline_vp) *
        data->timeline_scale;

      float viewport_t_offset = rig_ui_viewport_get_doc_x (data->timeline_vp);
      float viewport_y_offset = rig_ui_viewport_get_doc_y (data->timeline_vp);
      CoglPipeline *pipeline = cogl_pipeline_new (data->ctx->cogl_context);
      //NodeFloat *next;

      CoglPrimitive *prim;

      for (l = data->selected_transition->paths; l; l = l->next)
        {
          Path *path = l->data;
          GArray *points;
          GList *l;
          GList *next;
          float red, green, blue;

          if (path == NULL ||
              path->prop->object != item ||
              path->prop->spec->type != RIG_PROPERTY_TYPE_FLOAT)
            continue;

          if (strcmp (path->prop->spec->name, "x") == 0)
            red = 1.0, green = 0.0, blue = 0.0;
          else if (strcmp (path->prop->spec->name, "y") == 0)
            red = 0.0, green = 1.0, blue = 0.0;
          else if (strcmp (path->prop->spec->name, "z") == 0)
            red = 0.0, green = 0.0, blue = 1.0;
          else
            continue;

          points = g_array_new (FALSE, FALSE, sizeof (CoglVertexP2));

          for (l = path->nodes.head; l; l = next)
            {
              NodeFloat *f_node = l->data;
              CoglVertexP2 p;

              next = l->next;

              /* FIXME: This clipping wasn't working... */
#if 0
              /* Only draw the nodes within the current viewport */
              if (next)
                {
                  float max_t = viewport_t_offset + data->timeline_vp->width * viewport_t_scale;
                  if (next->t < viewport_t_offset)
                    continue;
                  if (node->t > max_t && next->t > max_t)
                    break;
                }
#endif

#define HANDLE_HALF_SIZE 4
              p.x = viewport_x + (f_node->t - viewport_t_offset) * viewport_t_scale;

              cogl_pipeline_set_color4f (pipeline, red, green, blue, 1);

              p.y = viewport_y + (f_node->value - viewport_y_offset) * viewport_y_scale;
#if 1
#if 1
              cogl_framebuffer_push_matrix (fb);
              cogl_framebuffer_translate (fb, p.x, p.y, 0);
              cogl_framebuffer_scale (fb, HANDLE_HALF_SIZE, HANDLE_HALF_SIZE, 0);
              cogl_framebuffer_draw_attributes (fb,
                                                pipeline,
                                                COGL_VERTICES_MODE_LINE_STRIP,
                                                1,
                                                data->circle_node_n_verts - 1,
                                                &data->circle_node_attribute,
                                                1);
              cogl_framebuffer_pop_matrix (fb);
#else
#if 0
              cogl_framebuffer_draw_rectangle (fb,
                                               pipeline,
                                               p.x - HANDLE_HALF_SIZE,
                                               p.y - HANDLE_HALF_SIZE,
                                               p.x + HANDLE_HALF_SIZE,
                                               p.y + HANDLE_HALF_SIZE);
#endif
#endif
#endif
              g_array_append_val (points, p);
            }

          prim = cogl_primitive_new_p2 (ctx, COGL_VERTICES_MODE_LINE_STRIP,
                                        points->len, (CoglVertexP2 *)points->data);
          draw_jittered_primitive4f (data, fb, prim, red, green, blue);
          cogl_object_unref (prim);

          g_array_free (points, TRUE);
        }

      cogl_object_unref (pipeline);

      {
        double progress;
        float progress_x;
        float progress_line[4];

        progress = rig_timeline_get_progress (data->timeline);

        progress_x = -viewport_t_offset * viewport_t_scale + data->timeline_width * data->timeline_scale * progress;
        progress_line[0] = progress_x;
        progress_line[1] = 0;
        progress_line[2] = progress_x;
        progress_line[3] = data->timeline_height;

        prim = cogl_primitive_new_p2 (ctx, COGL_VERTICES_MODE_LINE_STRIP, 2, (CoglVertexP2 *)progress_line);
        draw_jittered_primitive4f (data, fb, prim, 0, 1, 0);
        cogl_object_unref (prim);
      }
    }

  //cogl_framebuffer_pop_matrix (fb);

  rig_camera_end_frame (camera);
}

static RigTraverseVisitFlags
_rig_scenegraph_pre_paint_cb (RigObject *object,
                              int depth,
                              void *user_data)
{
  RigPaintContext *paint_ctx = user_data;
  RigCamera *camera = paint_ctx->camera;
  CoglFramebuffer *fb = rig_camera_get_framebuffer (camera);
  RigPaintableVTable *vtable =
    rig_object_get_vtable (object, RIG_INTERFACE_ID_PAINTABLE);

#if 0
  if (rig_object_get_type (object) == &rig_camera_type)
    {
      g_print ("%*sCamera = %p\n", depth, "", object);
      rig_camera_flush (RIG_CAMERA (object));
      return RIG_TRAVERSE_VISIT_CONTINUE;
    }
  else
#endif

  if (rig_object_get_type (object) == &rig_ui_viewport_type)
    {
      RigUIViewport *ui_viewport = RIG_UI_VIEWPORT (object);
#if 0
      g_print ("%*sPushing clip = %f %f\n",
               depth, "",
               rig_ui_viewport_get_width (ui_viewport),
               rig_ui_viewport_get_height (ui_viewport));
#endif
      cogl_framebuffer_push_rectangle_clip (fb,
                                            0, 0,
                                            rig_ui_viewport_get_width (ui_viewport),
                                            rig_ui_viewport_get_height (ui_viewport));
    }

  if (rig_object_is (object, RIG_INTERFACE_ID_TRANSFORMABLE))
    {
      //g_print ("%*sTransformable = %p\n", depth, "", object);
      const CoglMatrix *matrix = rig_transformable_get_matrix (object);
      //cogl_debug_matrix_print (matrix);
      cogl_framebuffer_push_matrix (fb);
      cogl_framebuffer_transform (fb, matrix);
    }

  if (rig_object_is (object, RIG_INTERFACE_ID_PAINTABLE))
    vtable->paint (object, paint_ctx);

  /* XXX:
   * How can we maintain state between the pre and post stages?  Is it
   * ok to just "sub-class" the paint context and maintain a stack of
   * state that needs to be shared with the post paint code.
   */

  return RIG_TRAVERSE_VISIT_CONTINUE;
}

static RigTraverseVisitFlags
_rig_scenegraph_post_paint_cb (RigObject *object,
                               int depth,
                               void *user_data)
{
  RigPaintContext *paint_ctx = user_data;
  CoglFramebuffer *fb = rig_camera_get_framebuffer (paint_ctx->camera);

#if 0
  if (rig_object_get_type (object) == &rig_camera_type)
    {
      rig_camera_end_frame (RIG_CAMERA (object));
      return RIG_TRAVERSE_VISIT_CONTINUE;
    }
  else
#endif

  if (rig_object_get_type (object) == &rig_ui_viewport_type)
    {
      cogl_framebuffer_pop_clip (fb);
    }

  if (rig_object_is (object, RIG_INTERFACE_ID_TRANSFORMABLE))
    {
      cogl_framebuffer_pop_matrix (fb);
    }

  return RIG_TRAVERSE_VISIT_CONTINUE;
}

static CoglBool
test_paint (RigShell *shell, void *user_data)
{
  Data *data = user_data;
  CoglFramebuffer *fb = rig_camera_get_framebuffer (data->camera);
  TestPaintContext test_paint_ctx;
  RigPaintContext *paint_ctx = &test_paint_ctx._parent;

  cogl_framebuffer_clear4f (fb,
                            COGL_BUFFER_BIT_COLOR|COGL_BUFFER_BIT_DEPTH,
                            0.5, 0.5, 0.5, 1);

  test_paint_ctx.data = data;
  test_paint_ctx.shadow_pass = FALSE;

  paint_ctx->camera = data->camera;

  //g_print ("Data camera = %p\n", data->camera);
  rig_camera_flush (data->camera);
  rig_graphable_traverse (data->root,
                          RIG_TRAVERSE_DEPTH_FIRST,
                          _rig_scenegraph_pre_paint_cb,
                          _rig_scenegraph_post_paint_cb,
                          &test_paint_ctx);
  rig_camera_end_frame (data->camera);

  paint_main_area_camera (data->main_camera, &test_paint_ctx);

  paint_timeline_camera (data->timeline_camera, data);

#if 0
  paint_ctx->camera = data->main_camera;

  rig_graphable_traverse (data->main_camera,
                          RIG_TRAVERSE_DEPTH_FIRST,
                          _rig_scenegraph_pre_paint_cb,
                          _rig_scenegraph_post_paint_cb,
                          &test_paint_ctx);
#endif

  cogl_onscreen_swap_buffers (COGL_ONSCREEN (fb));

  return FALSE;
}

#if 0
static void
path_t_update_cb (RigProperty *property, void *user_data)
{
  Data *data = user_data;
  double elapsed = rig_timeline_get_elapsed (data->timeline);
  Node *n0, *n1;
  float x, y, z;

  path_find_control_points2 (data->path, elapsed, 1, &n0, &n1);

  lerp_node_pos (n0, n1, elapsed, &x, &y, &z);

  g_print ("e=%f: 0(%f,%f,%f) 1(%f,%f,%f) lerp=(%f, %f, %f)\n", elapsed,
           n0->point[0], n0->point[1], n0->point[2],
           n1->point[0], n1->point[1], n1->point[2],
           x, y, z);
}
#endif

static RigIntrospectableVTable _asset_introspectable_vtable = {
  rig_simple_introspectable_lookup_property,
  rig_simple_introspectable_foreach_property
};

static RigType _asset_type;

static void
_asset_type_init (void)
{
  rig_type_init (&_asset_type);
#if 0
  rig_type_add_interface (&_asset_type,
                          RIG_INTERFACE_ID_INTROSPECTABLE,
                          0, /* no implied properties */
                          &_asset_introspectable_vtable);
  rig_type_add_interface (&_asset_type,
                          RIG_INTERFACE_ID_SIMPLE_INTROSPECTABLE,
                          offsetof (Asset, introspectable),
                          NULL); /* no implied vtable */
#endif
}

static Asset *
asset_new_texture (Data *data,
                   uint32_t id,
                   const char *path)
{
  Asset *asset = g_slice_new (Asset);
  char *full_path;
  //CoglTexture *texture;
  GError *error = NULL;

  rig_object_init (&asset->_parent, &_asset_type);

  asset->data = data;

  asset->id = id;

  asset->type = ASSET_TYPE_TEXTURE;

#ifndef __ANDROID__
  full_path = g_strconcat (_rig_project_dir, path, NULL);
  asset->texture = rig_load_texture (data->ctx, full_path, &error);
  g_free (full_path);
#else
  asset->texture = rig_load_texture (data->ctx, path, &error);
#endif

  if (!asset->texture)
    {
      g_slice_free (Asset, asset);
      g_error ("Failed to load asset texture: %s", error->message);
      g_error_free (error);
      return NULL;
    }

  asset->path = g_strdup (path);

  //rig_simple_introspectable_init (asset);

  return asset;
}

static void
_asset_free (Asset *asset)
{
  if (asset->texture)
    cogl_object_unref (asset->texture);

  //rig_simple_introspectable_destroy (asset);

  g_slice_free (Asset, asset);
}

static void
update_transition_progress_cb (RigProperty *property, void *user_data)
{
  Data *data = user_data;
  double elapsed = rig_timeline_get_elapsed (data->timeline);
  Transition *transition = property->object;
  //GList *l;

  transition->progress = elapsed;
  rig_property_dirty (&data->ctx->property_ctx,
                      &transition->props[TRANSITION_PROP_PROGRESS]);

#if 0
  /* FIXME: This should just update transition->progress and there
   * should be separate bindings to update all of the other properties
   */

  for (l = transition->paths; l; l = l->next)
    {
      Path *path = l->data;
      if (!path)
        continue;
      path_lerp_property (path, elapsed);
    }

  rig_transform_init_identity (item->transform);
  rig_transform_translate (item->transform,
                           item->x,
                           item->y,
                           item->z);
  rig_transform_quaternion_rotate (item->transform, &item->rotation);
#endif


  //cogl_matrix_init_identity (&item->transform);
  //cogl_matrix_translate (&item->transform, x, y, z);

#if 0
  g_print ("ITEM: e=%f: 0(%f,%f,%f) 1(%f,%f,%f) lerp=(%f, %f, %f)\n", elapsed,
           n0->point[0], n0->point[1], n0->point[2],
           n1->point[0], n1->point[1], n1->point[2],
           x, y, z);
#endif
}

static void
unproject_window_coord (RigCamera *camera,
                        const CoglMatrix *modelview,
                        const CoglMatrix *inverse_modelview,
                        float object_coord_z,
                        float *x,
                        float *y)
{
  const CoglMatrix *projection = rig_camera_get_projection (camera);
  const CoglMatrix *inverse_projection =
    rig_camera_get_inverse_projection (camera);
  //float z;
  float ndc_x, ndc_y, ndc_z, ndc_w;
  float eye_x, eye_y, eye_z, eye_w;
  const float *viewport = rig_camera_get_viewport (camera);

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

static RigInputEventStatus
item_grab_input_cb (RigInputEvent *event,
                    void *user_data)

{
  Item *item = user_data;
  Data *data = item->data;

  //g_print ("Item grab event\n");

  if (rig_input_event_get_type (event) == RIG_INPUT_EVENT_TYPE_MOTION)
    {
      float x = rig_motion_event_get_x (event);
      float y = rig_motion_event_get_y (event);
      CoglMatrix device_transform;
      CoglMatrix inverse_transform;
      RigCamera *camera = rig_input_event_get_camera (event);
      //double elapsed;

      rig_graphable_get_modelview (data->device_transform,
                                   camera,
                                   &device_transform);

      if (!cogl_matrix_get_inverse (&device_transform, &inverse_transform))
        g_error ("Failed to get inverse transform");

      unproject_window_coord (camera,
                              &device_transform,
                              &inverse_transform,
                              0, /* z in item coordinates */
                              &x, &y);

      if (rig_motion_event_get_action (event) == RIG_MOTION_EVENT_ACTION_UP)
        {
          Transition *transition = data->selected_transition;
          float elapsed = rig_timeline_get_elapsed (data->timeline);
          Path *path_x = transition_find_path (transition, "x");
          Path *path_y = transition_find_path (transition, "y");
          //Path *path_z = transition_find_path (transition, "z");


          path_insert_float (path_x, elapsed, x);
          path_insert_float (path_y, elapsed, y);
          //path_insert_float (item->prop_paths[ITEM_PROP_Z], elapsed, 0);

          rig_shell_ungrab_input (data->ctx->shell,
                                  item_grab_input_cb,
                                  user_data);

          rig_shell_queue_redraw (data->ctx->shell);

          return RIG_INPUT_EVENT_STATUS_HANDLED;
        }
      else if (rig_motion_event_get_action (event) == RIG_MOTION_EVENT_ACTION_MOVE)
        {
#if 1
          cogl_matrix_init_identity (&item->transform);
          cogl_matrix_translate (&item->transform, x, y, 0);
#else
          rig_transform_init_identity (item->transform);
          rig_transform_translate (item->transform, x, y, 0);
#endif
          rig_shell_queue_redraw (data->ctx->shell);

          return RIG_INPUT_EVENT_STATUS_HANDLED;
        }
    }

  return RIG_INPUT_EVENT_STATUS_UNHANDLED;
}

typedef struct _ToolListState
{
  Data *data;
  RigObject *tool_list;
  int n_items;
} ToolListState;

static void
property_updated_cb (RigText *text,
                     void *user_data)
{
  RigProperty *prop = user_data;
  RigContext *ctx = rig_text_get_context (text);

  switch (prop->spec->type)
    {
      case RIG_PROPERTY_TYPE_FLOAT:
        {
          const char *str = rig_text_get_text (text);
          float value = strtod (str, NULL);
          rig_property_set_float (&ctx->property_ctx, prop, value);
          break;
        }
      default:
        g_warning ("FIXME: missing code to sync text entry back to property");
    }

  rig_shell_queue_redraw (ctx->shell);
}

static void
add_item_property (RigProperty *prop, void *user_data)
{
  ToolListState *state = user_data;
  Data *data = state->data;
  //RigText *text;

  RigTransform *transform =
    rig_transform_new (data->ctx, NULL);

  switch (prop->spec->type)
    {
    case RIG_PROPERTY_TYPE_BOOLEAN:
      {
        RigToggle *toggle = rig_toggle_new (data->ctx,
                                            prop->spec->name);
        rig_graphable_add_child (transform, toggle);
        break;
      }
    case RIG_PROPERTY_TYPE_FLOAT:
      {
        RigText *label = rig_text_new (data->ctx);
        RigTransform *text_transform = rig_transform_new (data->ctx, NULL);
        RigText *text = rig_text_new (data->ctx);
        char *float_value = g_strdup_printf ("%f", rig_property_get_float (prop));
        float width, height;

        if (prop->spec->nick)
          rig_text_set_text (label, prop->spec->nick);
        else
          rig_text_set_text (label, prop->spec->name);

        rig_graphable_add_child (text_transform, text);

        rig_sizable_get_size (label, &width, &height);
        rig_transform_translate (text_transform, width + 10, 0, 0);

        rig_text_set_text (text, float_value);
        g_free (float_value);
        rig_text_set_color_u32 (text, 0xffffffff);
        rig_text_set_selection_color_u32 (text, 0x009ccfff);
        rig_text_set_editable (text, TRUE);
        rig_text_set_single_line_mode (text, TRUE);
        rig_text_set_activatable (text, TRUE);
        rig_text_set_activate_callback (text,
                                        property_updated_cb,
                                        prop);

        rig_graphable_add_child (transform, label);
        rig_graphable_add_child (transform, text_transform);
        break;
      }
    default:
      {
        RigText *label = rig_text_new (data->ctx);

        if (prop->spec->nick)
          rig_text_set_text (label, prop->spec->nick);
        else
          rig_text_set_text (label, prop->spec->name);

        rig_graphable_add_child (transform, label);
      }
    }

  rig_transform_translate (transform, 0, state->n_items++ * 50, 0);

  rig_graphable_add_child (state->tool_list, transform);
  g_print ("prop name = %s\n", prop->spec->name);
}

static void
update_tool_list (Data *data)
{
  //GList *l;
  //int i = 0;
  RigObject *doc_node;
  ToolListState state;

  if (data->tool_list)
    rig_graphable_remove_child (data->tool_list);

  data->tool_list = rig_graph_new (data->ctx, NULL);

  doc_node = rig_ui_viewport_get_doc_node (data->tool_vp);
  rig_graphable_add_child (doc_node, data->tool_list);
  rig_ref_countable_unref (data->tool_list);

  state.data = data;
  state.n_items = 0;
  state.tool_list = data->tool_list;

  rig_introspectable_foreach_property (data->selected_item, add_item_property, &state);
}


static RigInputEventStatus
item_input_region_cb (RigInputRegion *region,
                      RigInputEvent *event,
                      void *user_data)
{
  Item *item = user_data;
  Data *data = item->data;

  //g_print ("Item input\n");

  if (rig_input_event_get_type (event) == RIG_INPUT_EVENT_TYPE_MOTION)
    {
      if (rig_motion_event_get_action (event) == RIG_MOTION_EVENT_ACTION_DOWN)
        {
          data->selected_item = item;
          update_tool_list (data);
          rig_shell_grab_input (data->ctx->shell,
                                rig_input_event_get_camera (event),
                                item_grab_input_cb, item);
          return RIG_INPUT_EVENT_STATUS_HANDLED;
        }
    }

  return RIG_INPUT_EVENT_STATUS_UNHANDLED;
}

static RigIntrospectableVTable _item_introspectable_vtable = {
  rig_simple_introspectable_lookup_property,
  rig_simple_introspectable_foreach_property
};

static RigGraphableVTable _item_graphable_vtable = {
  NULL, /* child remove */
  NULL, /* child add */
  NULL /* parent changed */
};

static void
item_free (Item *item)
{
  cogl_object_unref (item->pipeline);
  //rig_ref_countable_unref (item->transform);

  rig_simple_introspectable_destroy (item);

  g_slice_free (Item, item);
}

RigRefCountableVTable _item_ref_countable_vtable = {
  rig_ref_countable_simple_ref,
  rig_ref_countable_simple_unref,
  item_free
};

static void
item_paint (RigObject *object,
            RigPaintContext *paint_ctx)
{
  Item *item = object;
#if 0
  RigCamera *camera = paint_ctx->camera;
  CoglFramebuffer *fb = rig_camera_get_framebuffer (camera);

  cogl_framebuffer_draw_rectangle (fb,
                                   item->pipeline,
                                   item->x0,
                                   item->y0,
                                   item->x1,
                                   item->y1);
#endif

  rig_paintable_paint (item->diamond_slice, paint_ctx);
}

static RigPaintableVTable _item_paintable_vtable = {
  item_paint
};

static const CoglMatrix *
item_get_matrix (Item *item)
{
  return &item->transform;
}

static RigTransformableVTable _item_transformable_vtable = {
  item_get_matrix
};

static RigType _item_type;

static void
_item_type_init (void)
{
  rig_type_init (&_item_type);
  rig_type_add_interface (&_item_type,
                          RIG_INTERFACE_ID_REF_COUNTABLE,
                          offsetof (Item, ref_count),
                          &_item_ref_countable_vtable);
  rig_type_add_interface (&_item_type,
                          RIG_INTERFACE_ID_INTROSPECTABLE,
                          0, /* no implied properties */
                          &_item_introspectable_vtable);
  rig_type_add_interface (&_item_type,
                          RIG_INTERFACE_ID_SIMPLE_INTROSPECTABLE,
                          offsetof (Item, introspectable),
                          NULL); /* no implied vtable */
  rig_type_add_interface (&_item_type,
                          RIG_INTERFACE_ID_GRAPHABLE,
                          offsetof (Item, graphable),
                          &_item_graphable_vtable);
  rig_type_add_interface (&_item_type,
                          RIG_INTERFACE_ID_PAINTABLE,
                          offsetof (Item, paintable),
                          &_item_paintable_vtable);
  rig_type_add_interface (&_item_type,
                           RIG_INTERFACE_ID_TRANSFORMABLE,
                           0,
                           &_item_transformable_vtable);
}

static void
update_item_transform_cb (RigProperty *transform_prop, void *user_data)
{
  Item *item = transform_prop->object;

  /* FIXME: Add rig_transform_init_translate() */

  cogl_matrix_init_translation (&item->transform,
                                item->x,
                                item->y,
                                item->z);
  cogl_matrix_rotate_quaternion (&item->transform,
                                 &item->rotation);
#if 0
  rig_transform_init_identity (item->transform);
  rig_transform_translate (item->transform,
                           item->x,
                           item->y,
                           item->z);
  rig_transform_quaternion_rotate (item->transform, &item->rotation);
#endif
}

Item *
item_new (Data *data,
          Asset *texture_asset,
          uint32_t id)
{
  CoglTexture *texture;
  Item *item;
  float width, height;
  float half_diamond_size;

  g_return_val_if_fail (texture_asset->type == ASSET_TYPE_TEXTURE, NULL);

  item = g_slice_new0 (Item);

  item->ref_count = 1;

  item->id = id;

  rig_object_init (&item->_parent, &_item_type);

  rig_graphable_init (item);
  rig_paintable_init (item);

  rig_simple_introspectable_init (item, _item_prop_specs, item->props);

  texture = texture_asset->texture;

  item->diamond_size = 400;
  half_diamond_size = item->diamond_size / 2.0;

  item->diamond_slice = diamond_slice_new (data,
                                           texture,
                                           item->diamond_size);

  item->texture_asset = texture_asset;

  item->data = data;

  item->x = item->y = item->z = 0;
  cogl_quaternion_init (&item->rotation, 0, 0, 1, 0);
  cogl_matrix_init_identity (&item->transform);

  item->pipeline = cogl_pipeline_new (data->ctx->cogl_context);
  cogl_pipeline_set_layer_texture (item->pipeline, 0, texture);

  //item->x0 = 0;
  //item->y0 = 0;
  width = cogl_texture_get_width (texture);
  height = cogl_texture_get_height (texture);
  //item->x1 = width;
  //item->y1 = height;

  item->input_transform =
    rig_transform_new (data->ctx,
                       item->input_region =
                       rig_input_region_new_rectangle (0, 0,
                                                       item->diamond_size,
                                                       item->diamond_size,
                                                       item_input_region_cb,
                                                       item),
                       NULL);

  rig_transform_rotate (item->input_transform, 45, 0, 0, 1);
  rig_transform_translate (item->input_transform, - half_diamond_size, - half_diamond_size, 0);
  //rig_transform_translate (item->input_transform, half_diamond_size, half_diamond_size, 0);

  rig_graphable_add_child (item, item->input_transform);
  //rig_input_region_set_graphable (item->input_region, item->input_transform);

  //cogl_matrix_init_identity (&item->transform_matrix);
  //rig_input_region_set_transform (item->input_region, &item->transform_matrix);
  //rig_camera_add_input_region (data->main_camera, item->input_region);

  rig_property_set_binding (&item->props[ITEM_PROP_TRANSFORM],
                            update_item_transform_cb,
                            NULL,
                            &item->props[ITEM_PROP_X],
                            &item->props[ITEM_PROP_Y],
                            &item->props[ITEM_PROP_Z],
                            &item->props[ITEM_PROP_ROTATION],
                            NULL);
  return item;
}

static RigIntrospectableVTable _transition_introspectable_vtable = {
  rig_simple_introspectable_lookup_property,
  rig_simple_introspectable_foreach_property
};

static RigType _transition_type;

static void
_transition_type_init (void)
{
  rig_type_init (&_transition_type);
  rig_type_add_interface (&_transition_type,
                          RIG_INTERFACE_ID_INTROSPECTABLE,
                          0, /* no implied properties */
                          &_transition_introspectable_vtable);
  rig_type_add_interface (&_transition_type,
                          RIG_INTERFACE_ID_SIMPLE_INTROSPECTABLE,
                          offsetof (Transition, introspectable),
                          NULL); /* no implied vtable */
}

Transition *
transition_new (Data *data,
                uint32_t id)
{
  //GError *error = NULL;
  Transition *transition;

  transition = g_slice_new0 (Transition);

  transition->id = id;

  rig_object_init (&transition->_parent, &_transition_type);

  rig_simple_introspectable_init (transition, _transition_prop_specs, transition->props);

  transition->data = data;

  transition->progress = 0;
  transition->paths = NULL;

#if 0
  rig_property_init (&transition->props[ITEM_PROP_PROGRESS],
                     NULL, /* dummy property: no spec */
                     transition);
#endif

  rig_property_set_binding (&transition->props[TRANSITION_PROP_PROGRESS],
                            update_transition_progress_cb,
                            data,
                            data->timeline_elapsed,
                            NULL);

  return transition;
}

static void
transition_free (Transition *transition)
{
  GList *l;

  rig_simple_introspectable_destroy (transition);

  for (l = transition->paths; l; l = l->next)
    {
      Path *path = l->data;
      path_free (path);
    }

  g_slice_free (Transition, transition);
}

void
transition_add_path (Transition *transition,
                     Path *path)
{
  transition->paths = g_list_prepend (transition->paths, path);
}

static Path *
transition_find_path (Transition *transition,
                      const char *property)
{
  GList *l;

  for (l = transition->paths; l; l = l->next)
    {
      Path *path = l->data;

      if (strcmp (path->prop->spec->name, property) == 0)
        return path;
    }

  return NULL;
}
#if 0
static void
update_slider_pos_cb (RigProperty *property,
                      void *user_data)
{
  Data *data = user_data;
  double progress = rig_timeline_get_progress (data->timeline);

  //g_print ("update_slider_pos_cb %f\n", progress);

  rig_slider_set_progress (data->slider, progress);
}

static void
update_timeline_progress_cb (RigProperty *property,
                             void *user_data)
{
  Data *data = user_data;
  double progress = rig_property_get_float (data->slider_progress);

  //g_print ("update_timeline_progress_cb %f\n", progress);

  rig_timeline_set_progress (data->timeline, progress);
}
#endif

static RigInputEventStatus
timeline_grab_input_cb (RigInputEvent *event, void *user_data)
{
  Data *data = user_data;

  if (rig_input_event_get_type (event) != RIG_INPUT_EVENT_TYPE_MOTION)
    return RIG_INPUT_EVENT_STATUS_UNHANDLED;

  if (rig_motion_event_get_action (event) == RIG_MOTION_EVENT_ACTION_MOVE)
    {
      RigButtonState state = rig_motion_event_get_button_state (event);
      float x = rig_motion_event_get_x (event);
      float y = rig_motion_event_get_y (event);

      if (state & RIG_BUTTON_STATE_1)
        {
          RigCamera *camera = data->timeline_camera;
          const CoglMatrix *view = rig_camera_get_view_transform (camera);
          CoglMatrix inverse_view;
          float progress;

          if (!cogl_matrix_get_inverse (view, &inverse_view))
            g_error ("Failed to get inverse transform");

          unproject_window_coord (camera,
                                  view,
                                  &inverse_view,
                                  0, /* z in item coordinates */
                                  &x, &y);

          progress = x / data->timeline_width;
          //g_print ("x = %f, y = %f progress=%f\n", x, y, progress);

          rig_timeline_set_progress (data->timeline, progress);
          rig_shell_queue_redraw (data->ctx->shell);

          return RIG_INPUT_EVENT_STATUS_HANDLED;
        }
      else if (state & RIG_BUTTON_STATE_2)
        {
          float dx = data->grab_x - x;
          float dy = data->grab_y - y;
          float t_scale =
            rig_ui_viewport_get_doc_scale_x (data->timeline_vp) *
            data->timeline_scale;

          float y_scale =
            rig_ui_viewport_get_doc_scale_y (data->timeline_vp) *
            data->timeline_scale;

          float inv_t_scale = 1.0 / t_scale;
          float inv_y_scale = 1.0 / y_scale;


          rig_ui_viewport_set_doc_x (data->timeline_vp,
                                     data->grab_timeline_vp_t + (dx * inv_t_scale));
          rig_ui_viewport_set_doc_y (data->timeline_vp,
                                     data->grab_timeline_vp_y + (dy * inv_y_scale));

          rig_shell_queue_redraw (data->ctx->shell);
        }
    }
  else if (rig_motion_event_get_action (event) == RIG_MOTION_EVENT_ACTION_UP)
    {
      rig_shell_ungrab_input (data->ctx->shell,
                              timeline_grab_input_cb,
                              user_data);
      return RIG_INPUT_EVENT_STATUS_HANDLED;
    }
  return RIG_INPUT_EVENT_STATUS_UNHANDLED;
}

static RigInputEventStatus
timeline_input_cb (RigInputEvent *event,
                   void *user_data)
{
  Data *data = user_data;

  if (rig_input_event_get_type (event) == RIG_INPUT_EVENT_TYPE_MOTION)
    {
      data->key_focus_callback = timeline_input_cb;

      switch (rig_motion_event_get_action (event))
        {
        case RIG_MOTION_EVENT_ACTION_DOWN:
          {
            data->grab_x = rig_motion_event_get_x (event);
            data->grab_y = rig_motion_event_get_y (event);
            data->grab_timeline_vp_t = rig_ui_viewport_get_doc_x (data->timeline_vp);
            data->grab_timeline_vp_y = rig_ui_viewport_get_doc_y (data->timeline_vp);
            /* TODO: Add rig_shell_implicit_grab_input() that handles releasing
             * the grab for you */
            g_print ("timeline input grab\n");
            rig_shell_grab_input (data->ctx->shell,
                                  rig_input_event_get_camera (event),
                                  timeline_grab_input_cb, data);
            return RIG_INPUT_EVENT_STATUS_HANDLED;
          }
	default:
	  break;
        }
    }
  else if (rig_input_event_get_type (event) == RIG_INPUT_EVENT_TYPE_KEY &&
           rig_key_event_get_action (event) == RIG_KEY_EVENT_ACTION_UP)
    {
      switch (rig_key_event_get_keysym (event))
        {
        case RIG_KEY_equal:
          data->timeline_scale += 0.2;
          rig_shell_queue_redraw (data->ctx->shell);
          break;
        case RIG_KEY_minus:
          data->timeline_scale -= 0.2;
          rig_shell_queue_redraw (data->ctx->shell);
          break;
        case RIG_KEY_Home:
          data->timeline_scale = 1;
          rig_shell_queue_redraw (data->ctx->shell);
        }
      g_print ("Key press in timeline area\n");
    }

  return RIG_INPUT_EVENT_STATUS_UNHANDLED;
}

static RigInputEventStatus
timeline_region_input_cb (RigInputRegion *region,
                          RigInputEvent *event,
                          void *user_data)
{
  return timeline_input_cb (event, user_data);
}

static CoglPrimitive *
create_line_primitive (float a[3], float b[3])
{
  CoglVertexP3 data[2];
  CoglAttributeBuffer *attribute_buffer;
  CoglAttribute *attributes[1];
  CoglPrimitive *primitive;

  data[0].x = a[0];
  data[0].y = a[1];
  data[0].z = a[2];
  data[1].x = b[0];
  data[1].y = b[1];
  data[1].z = b[2];

  attribute_buffer = cogl_attribute_buffer_new (rig_cogl_context,
                                                2 * sizeof (CoglVertexP3),
                                                data);

  attributes[0] = cogl_attribute_new (attribute_buffer,
                                      "cogl_position_in",
                                      sizeof (CoglVertexP3),
                                      offsetof (CoglVertexP3, x),
                                      3,
                                      COGL_ATTRIBUTE_TYPE_FLOAT);

  primitive = cogl_primitive_new_with_attributes (COGL_VERTICES_MODE_LINES,
                                                  2, attributes, 1);

  cogl_object_unref (attribute_buffer);
  cogl_object_unref (attributes[0]);

  return primitive;
}

static void
transform_ray (CoglMatrix *transform,
               bool        inverse_transform,
               float       ray_origin[3],
               float       ray_direction[3])
{
  CoglMatrix inverse, normal_matrix, *m;

  m = transform;
  if (inverse_transform)
    {
      cogl_matrix_get_inverse (transform, &inverse);
      m = &inverse;
    }

  cogl_matrix_transform_points (m,
                                3, /* num components for input */
                                sizeof (float) * 3, /* input stride */
                                ray_origin,
                                sizeof (float) * 3, /* output stride */
                                ray_origin,
                                1 /* n_points */);

  cogl_matrix_get_inverse (m, &normal_matrix);
  cogl_matrix_transpose (&normal_matrix);

  rig_util_transform_normal (&normal_matrix,
                             &ray_direction[0],
                             &ray_direction[1],
                             &ray_direction[2]);
}

static CoglPrimitive *
create_picking_ray (Data            *data,
                    CoglFramebuffer *fb,
                    float            ray_position[3],
                    float            ray_direction[3],
                    float            length)
{
  CoglPrimitive *line;
  float points[6];

  points[0] = ray_position[0];
  points[1] = ray_position[1];
  points[2] = ray_position[2];
  points[3] = ray_position[0] + length * ray_direction[0];
  points[4] = ray_position[1] + length * ray_direction[1];
  points[5] = ray_position[2] + length * ray_direction[2];

  line = create_line_primitive (points, points + 3);

  return line;
}

static RigEntity *
pick (Data *data,
      float ray_origin[3],
      float ray_direction[3])
{
  RigEntity *entity, *selected_entity = NULL;
  RigComponent *mesh;
  uint8_t *vertex_data;
  int n_vertices;
  size_t stride;
  int index;
  float distance, min_distance = G_MAXFLOAT;
  bool hit;
  float transformed_ray_origin[3];
  float transformed_ray_direction[3];
  GList *l;

  for (l = data->pickables; l; l = l->next)
    {
      entity = l->data;

      /* transform the ray into the model space */
      memcpy (transformed_ray_origin, ray_origin, 3 * sizeof (float));
      memcpy (transformed_ray_direction, ray_direction, 3 * sizeof (float));
      transform_ray (rig_entity_get_transform (entity),
                     TRUE, /* inverse of the transform */
                     transformed_ray_origin,
                     transformed_ray_direction);

      /* intersect the transformed ray with the mesh data */
      mesh = rig_entity_get_component (entity, RIG_COMPONENT_TYPE_GEOMETRY);
      if (!mesh || rig_object_get_type (mesh) != &rig_mesh_renderer_type)
        continue;

      vertex_data =
        rig_mesh_renderer_get_vertex_data (RIG_MESH_RENDERER (mesh), &stride);
      if (!vertex_data)
        continue;

      n_vertices = rig_mesh_renderer_get_n_vertices (RIG_MESH_RENDERER (mesh));

      hit = rig_util_intersect_mesh (vertex_data,
                                     n_vertices,
                                     stride,
                                     transformed_ray_origin,
                                     transformed_ray_direction,
                                     &index,
                                     &distance);

      /* to compare intersection distances we need to re-transform it back
       * to the world space, */
      cogl_vector3_normalize (transformed_ray_direction);
      transformed_ray_direction[0] *= distance;
      transformed_ray_direction[1] *= distance;
      transformed_ray_direction[2] *= distance;

      rig_util_transform_normal (rig_entity_get_transform (entity),
                                 &transformed_ray_direction[0],
                                 &transformed_ray_direction[1],
                                 &transformed_ray_direction[2]);
      distance = cogl_vector3_magnitude (transformed_ray_direction);

      if (hit && distance < min_distance)
        {
          min_distance = distance;
          selected_entity = entity;
        }
    }

  if (selected_entity)
    {
      g_message ("Hit entity, triangle #%d, distance %.2f",
                 index, distance);
    }

  return selected_entity;
}

static void
update_camera_position (Data *data)
{
#if 0
#if 0
  /* Calculate where the origin currently is from the camera's
   * point of view. Then we can fixup the camera's position
   * so this matches the real position of the origin. */
  float relative_origin[3] = { 0, 0, -data->main_camera_z };


  rig_entity_get_transformed_position (data->main_camera_rig,
                                       relative_origin);
#else
  float relative_origin[3] = { data->main_width / 2, data->main_height / 2, 0 };
  RigCamera *camera_component =
    rig_entity_get_component (data->main_camera, RIG_COMPONENT_TYPE_CAMERA);
  const CoglMatrix *view;
  float w = 1;

  camera_update_view (data, data->main_camera, FALSE);
  view = rig_camera_get_view_transform (data->main_camera);

  cogl_matrix_transform_point (view,
                               &relative_origin[0],
                               &relative_origin[1],
                               &relative_origin[2],
                               &w);
#endif

  rig_entity_translate (data->main_camera_rig0,
                        data->origin[0] - relative_origin[0],
                        data->origin[1] - relative_origin[1],
                        data->origin[2] - relative_origin[2]);

#endif
  float pos[3] = {0, 0, 0};
  rig_entity_set_position (data->main_camera_rig0, pos);
  rig_entity_translate (data->main_camera_rig0,
                        data->origin[0],
                        data->origin[1],
                        data->origin[2]);

  rig_entity_set_position (data->main_camera_rig2, pos);
  rig_entity_translate (data->main_camera_rig2,
                        -data->origin[0],
                        -data->origin[1],
                        -data->origin[2]);
  rig_shell_queue_redraw (data->ctx->shell);
}

static void
print_quaternion (const CoglQuaternion *q,
                  const char *label)
{
  float angle = cogl_quaternion_get_rotation_angle (q);
  float axis[3];
  cogl_quaternion_get_rotation_axis (q, axis);
  g_print ("%s: [%f (%f, %f, %f)]\n", label, angle, axis[0], axis[1], axis[2]);
}

static RigInputEventStatus
main_input_cb (RigInputEvent *event,
               void *user_data)
{
  Data *data = user_data;

  g_print ("Main Input Callback\n");

  if (rig_input_event_get_type (event) == RIG_INPUT_EVENT_TYPE_MOTION)
    {
      RigMotionEventAction action = rig_motion_event_get_action (event);
      RigModifierState modifiers = rig_motion_event_get_modifier_state (event);
      float x = rig_motion_event_get_x (event);
      float y = rig_motion_event_get_y (event);
      RigButtonState state;

      if (rig_camera_transform_window_coordinate (data->main_camera_component, &x, &y))
        data->key_focus_callback = main_input_cb;

      state = rig_motion_event_get_button_state (event);

      if (action == RIG_MOTION_EVENT_ACTION_DOWN &&
          state == RIG_BUTTON_STATE_1)
        {
          /* pick */
          RigComponent *camera;
          float ray_position[3], ray_direction[3], screen_pos[2],
                z_far, z_near;
          const float *viewport;
          const CoglMatrix *inverse_projection;
          //CoglMatrix *camera_transform;
          const CoglMatrix *camera_view;
          CoglMatrix camera_transform;

          camera = rig_entity_get_component (data->main_camera,
                                             RIG_COMPONENT_TYPE_CAMERA);
          viewport = rig_camera_get_viewport (RIG_CAMERA (camera));
          z_near = rig_camera_get_near_plane (RIG_CAMERA (camera));
          z_far = rig_camera_get_far_plane (RIG_CAMERA (camera));
          inverse_projection =
            rig_camera_get_inverse_projection (RIG_CAMERA (camera));

#if 0
          camera_transform = rig_entity_get_transform (data->main_camera);
#else
          camera_view = rig_camera_get_view_transform (camera);
          cogl_matrix_get_inverse (camera_view, &camera_transform);
#endif

          screen_pos[0] = x;
          screen_pos[1] = y;

          rig_util_create_pick_ray (viewport,
                                    inverse_projection,
                                    &camera_transform,
                                    screen_pos,
                                    ray_position,
                                    ray_direction);

          if (data->debug_pick_ray)
            {
              float x1 = 0, y1 = 0, z1 = z_near, w1 = 1;
              float x2 = 0, y2 = 0, z2 = z_far, w2 = 1;
              float len;

              if (data->picking_ray)
                cogl_object_unref (data->picking_ray);

              /* FIXME: This is a hack, we should intersect the ray with
               * the far plane to decide how long the debug primitive
               * should be */
              cogl_matrix_transform_point (&camera_transform,
                                           &x1, &y1, &z1, &w1);
              cogl_matrix_transform_point (&camera_transform,
                                           &x2, &y2, &z2, &w2);
              len = z2 - z1;

              data->picking_ray = create_picking_ray (data,
                                                      rig_camera_get_framebuffer (camera),
                                                      ray_position,
                                                      ray_direction,
                                                      len);
            }

          data->selected_entity = pick (data, ray_position, ray_direction);
          rig_shell_queue_redraw (data->ctx->shell);
          if (data->selected_entity == NULL)
            rig_tool_update (data->tool, NULL);
        }
      else if (action == RIG_MOTION_EVENT_ACTION_DOWN &&
          state == RIG_BUTTON_STATE_2)
        {
          //data->saved_rotation = *rig_entity_get_rotation (data->main_camera);
          data->saved_rotation = *rig_entity_get_rotation (data->main_camera_rig1);

          cogl_quaternion_init_identity (&data->arcball.q_drag);

          //rig_arcball_mouse_down (&data->arcball, data->width - x, y);
          rig_arcball_mouse_down (&data->arcball, data->main_width - x, data->main_height - y);
          g_print ("Arcball init, mouse = (%d, %d)\n", (int)(data->width - x), (int)(data->height - y));

          print_quaternion (&data->saved_rotation, "Saved Quaternion");
          print_quaternion (&data->arcball.q_drag, "Arcball Initial Quaternion");
          //data->button_down = TRUE;

          data->grab_x = x;
          data->grab_y = y;
          memcpy (data->saved_origin, data->origin, sizeof (data->origin));

          return RIG_INPUT_EVENT_STATUS_HANDLED;
        }
      else if (action == RIG_MOTION_EVENT_ACTION_MOVE &&
               state == RIG_BUTTON_STATE_2 &&
               modifiers & RIG_MODIFIER_SHIFT_ON)
        {
          float origin[3] = {0, 0, 0};
          float unit_x[3] = {1, 0, 0};
          float unit_y[3] = {0, 1, 0};
          float x_vec[3];
          float y_vec[3];
          float dx;
          float dy;
          float translation[3];

          rig_entity_get_transformed_position (data->main_camera,
                                               origin);
          rig_entity_get_transformed_position (data->main_camera,
                                               unit_x);
          rig_entity_get_transformed_position (data->main_camera,
                                               unit_y);

          x_vec[0] = origin[0] - unit_x[0];
          x_vec[1] = origin[1] - unit_x[1];
          x_vec[2] = origin[2] - unit_x[2];

          y_vec[0] = origin[0] - unit_y[0];
          y_vec[1] = origin[1] - unit_y[1];
          y_vec[2] = origin[2] - unit_y[2];

          dx = (x - data->grab_x) * (data->main_camera_z / 100.0f);
          dy = -(y - data->grab_y) * (data->main_camera_z / 100.0f);

          translation[0] = dx * x_vec[0];
          translation[1] = dx * x_vec[1];
          translation[2] = dx * x_vec[2];

          translation[0] += dy * y_vec[0];
          translation[1] += dy * y_vec[1];
          translation[2] += dy * y_vec[2];

          data->origin[0] = data->saved_origin[0] + translation[0];
          data->origin[1] = data->saved_origin[1] + translation[1];
          data->origin[2] = data->saved_origin[2] + translation[2];

          update_camera_position (data);

          g_print ("Translate %f %f dx=%f, dy=%f\n",
                   x - data->grab_x,
                   y - data->grab_y,
                   dx, dy);

          return RIG_INPUT_EVENT_STATUS_HANDLED;
        }
      else if (action == RIG_MOTION_EVENT_ACTION_MOVE &&
               state == RIG_BUTTON_STATE_2)
        {
          CoglQuaternion new_rotation;

          //if (!data->button_down)
          //  break;

          //rig_arcball_mouse_motion (&data->arcball, data->width - x, y);
          rig_arcball_mouse_motion (&data->arcball, data->main_width - x, data->main_height - y);
          g_print ("Arcball motion, center=%f,%f mouse = (%f, %f)\n",
                   data->arcball.center[0],
                   data->arcball.center[1],
                   x, y);

          cogl_quaternion_multiply (&new_rotation,
                                    &data->saved_rotation,
                                    &data->arcball.q_drag);

          //rig_entity_set_rotation (data->main_camera, &new_rotation);
          rig_entity_set_rotation (data->main_camera_rig1, &new_rotation);

          print_quaternion (&new_rotation, "New Rotation");

          print_quaternion (&data->arcball.q_drag, "Arcball Quaternion");

          g_print ("rig entity set rotation\n");
          /* XXX: The remaining problem is calculating the new
           * position for the camera!
           *
           * If we transform the point (0, 0, camera_z) by the
           * camera's transform we can find where the origin is
           * relative to the camera, and then find out how far that
           * point is from the true origin so we know how to
           * translate the camera.
           */
          update_camera_position (data);

          //rig_shell_queue_redraw (data->ctx->shell);

          return RIG_INPUT_EVENT_STATUS_HANDLED;
        }
      else if (action == RIG_MOTION_EVENT_ACTION_MOVE &&
               state == RIG_BUTTON_STATE_2 &&
               modifiers & RIG_MODIFIER_SHIFT_ON)
        {
          g_print ("Translate\n");
          return RIG_INPUT_EVENT_STATUS_HANDLED;
        }

#if 0
      switch (rig_motion_event_get_action (event))
        {
        case RIG_MOTION_EVENT_ACTION_DOWN:
          /* TODO: Add rig_shell_implicit_grab_input() that handles releasing
           * the grab for you */
          rig_shell_grab_input (data->ctx->shell, timeline_grab_input_cb, data);
          return RIG_INPUT_EVENT_STATUS_HANDLED;
        }
#endif
    }
  else if (rig_input_event_get_type (event) == RIG_INPUT_EVENT_TYPE_KEY &&
           rig_key_event_get_action (event) == RIG_KEY_EVENT_ACTION_UP)
    {
      switch (rig_key_event_get_keysym (event))
        {
        case RIG_KEY_s:
          save (data);
          break;
        case RIG_KEY_minus:
          if (data->main_camera_z)
            data->main_camera_z *= 1.2f;
          else
            data->main_camera_z = 0.1;

          update_camera_position (data);
          break;
        case RIG_KEY_equal:
          data->main_camera_z *= 0.8;
          update_camera_position (data);
          break;
        }
    }

  return RIG_INPUT_EVENT_STATUS_UNHANDLED;
}

static RigInputEventStatus
main_input_region_cb (RigInputRegion *region,
                      RigInputEvent *event,
                      void *user_data)
{
  return main_input_cb (event, user_data);
}

static void
allocate (Data *data)
{
  float screen_aspect;
  float main_aspect;
  float device_scale;
  float vp_width;
  float vp_height;

  data->top_bar_height = 30;
  //data->top_bar_height = 0;
  data->left_bar_width = data->width * 0.2;
  //data->left_bar_width = 200;
  //data->left_bar_width = 0;
  data->right_bar_width = data->width * 0.2;
  //data->right_bar_width = 200;
  data->bottom_bar_height = data->height * 0.2;
  data->grab_margin = 5;
  data->main_width = data->width - data->left_bar_width - data->right_bar_width;
  data->main_height = data->height - data->top_bar_height - data->bottom_bar_height;


  /* Update the window camera */
  rig_camera_set_projection_mode (data->camera, RIG_PROJECTION_ORTHOGRAPHIC);
  rig_camera_set_orthographic_coordinates (data->camera,
                                           0, 0, data->width, data->height);
  rig_camera_set_near_plane (data->camera, -1);
  rig_camera_set_far_plane (data->camera, 100);

  rig_camera_set_viewport (data->camera, 0, 0, data->width, data->height);

  screen_aspect = DEVICE_WIDTH / DEVICE_HEIGHT;
  main_aspect = data->main_width / data->main_height;

  rig_transform_init_identity (data->screen_area_transform);

  if (screen_aspect < main_aspect) /* screen is slimmer and taller than the main area */
    {
      data->screen_area_height = data->main_height;
      data->screen_area_width = data->screen_area_height * screen_aspect;

      rig_transform_translate (data->screen_area_transform,
                               (data->main_width / 2.0) - (data->screen_area_width / 2.0),
                               0, 0);
    }
  else
    {
      data->screen_area_width = data->main_width;
      data->screen_area_height = data->screen_area_width / screen_aspect;

      rig_transform_translate (data->screen_area_transform,
                               0,
                               (data->main_height / 2.0) - (data->screen_area_height / 2.0),
                               0);
    }

  rig_transform_init_identity (data->device_transform);

  /* NB: We know the screen area matches the device aspect ratio so we can use
   * a uniform scale here... */
  device_scale = data->screen_area_width / DEVICE_WIDTH;
  rig_transform_scale (data->device_transform,
                       device_scale,
                       device_scale,
                       device_scale);

  /* Setup projection for main content view */
  {
    float fovy = 10; /* y-axis field of view */
    float aspect = (float)data->main_width/(float)data->main_height;
    float z_near = 10; /* distance to near clipping plane */
    float z_far = 100; /* distance to far clipping plane */
#if 1
    fovy = 60;
    z_near = 1.1;
    z_far = 100;
#endif
    CoglMatrix projection;

    data->z_2d = 30; /* position to 2d plane */

    cogl_matrix_init_identity (&data->main_view);
    cogl_matrix_view_2d_in_perspective (&data->main_view,
                                        fovy, aspect, z_near, data->z_2d,
                                        data->main_width,
                                        data->main_height);
#if 0
    rig_camera_set_view_transform (data->main_camera_component, &data->main_view);
#endif

    rig_camera_set_projection_mode (data->main_camera_component,
                                    RIG_PROJECTION_PERSPECTIVE);
    rig_camera_set_field_of_view (data->main_camera_component, fovy);
    rig_camera_set_near_plane (data->main_camera_component, z_near);
    rig_camera_set_far_plane (data->main_camera_component, z_far);

#if 0
    cogl_matrix_init_identity (&data->main_view);
    rig_camera_set_projection_mode (data->main_camera_component,
                                    RIG_PROJECTION_ORTHOGRAPHIC);
    rig_camera_set_orthographic_coordinates (data->main_camera_component,
                                             0, 0, data->main_width, data->main_height);
    rig_camera_set_near_plane (data->main_camera_component, -1);
    rig_camera_set_far_plane (data->main_camera_component, 100);
#endif

    rig_camera_set_viewport (data->main_camera_component,
                             data->left_bar_width,
                             data->top_bar_height,
                             data->main_width,
                             data->main_height);

    rig_input_region_set_rectangle (data->main_input_region,
                                    data->left_bar_width,
                                    data->top_bar_height,
                                    data->left_bar_width + data->main_width,
                                    data->top_bar_height + data->main_height);
  }

#if 0
  data->origin[0] = data->main_width / 2;
  data->origin[1] = data->main_height / 2;
  data->origin[2] = 0;
  //data->origin[2] = data->z_2d;
#endif

  /* Setup projection for the timeline view */
  {
    data->timeline_width = data->width - data->right_bar_width;
    data->timeline_height = data->bottom_bar_height;

    rig_camera_set_projection_mode (data->timeline_camera, RIG_PROJECTION_ORTHOGRAPHIC);
    rig_camera_set_orthographic_coordinates (data->timeline_camera,
                                             0, 0,
                                             data->timeline_width,
                                             data->timeline_height);
    rig_camera_set_near_plane (data->timeline_camera, -1);
    rig_camera_set_far_plane (data->timeline_camera, 100);
    rig_camera_set_background_color4f (data->timeline_camera, 1, 0, 0, 1);

    rig_camera_set_viewport (data->timeline_camera,
                             0,
                             data->height - data->bottom_bar_height,
                             data->timeline_width,
                             data->timeline_height);

    rig_input_region_set_rectangle (data->timeline_input_region,
                                    0, 0, data->timeline_width, data->timeline_height);
  }

  vp_width = data->width - data->bottom_bar_height;
  rig_ui_viewport_set_width (data->timeline_vp,
                             vp_width);
  vp_height = data->bottom_bar_height;
  rig_ui_viewport_set_height (data->timeline_vp,
                              vp_height);
  rig_ui_viewport_set_doc_scale_x (data->timeline_vp,
                                   (vp_width / data->timeline_len));
  rig_ui_viewport_set_doc_scale_y (data->timeline_vp,
                                   (vp_height / DEVICE_HEIGHT));


#if 0
  {
    CoglMatrix input_transform;
    cogl_matrix_init_identity (&input_transform);
    cogl_matrix_translate (&input_transform, 0, -data->top_bar_height, 0);

    //rig_camera_set_input_transform (data->main_camera_component, &input_transform);
  }
#endif

  rig_rectangle_set_width (data->top_bar_rect, data->width);
  rig_rectangle_set_height (data->top_bar_rect, data->top_bar_height - data->grab_margin);

  {
    float left_bar_height = data->height - data->top_bar_height - data->bottom_bar_height;

    rig_rectangle_set_width (data->left_bar_rect, data->left_bar_width - data->grab_margin);
    rig_rectangle_set_height (data->left_bar_rect,
                              data->height - data->top_bar_height - data->bottom_bar_height);

    rig_transform_init_identity (data->left_bar_transform);
    rig_transform_translate (data->left_bar_transform,
                             0, data->top_bar_height, 0);

    rig_ui_viewport_set_width (data->assets_vp, data->left_bar_width);
    rig_ui_viewport_set_height (data->assets_vp, left_bar_height);
  }

  {
    float right_bar_height = data->height - data->top_bar_height;

    rig_rectangle_set_width (data->right_bar_rect, data->right_bar_width - data->grab_margin);
    rig_rectangle_set_height (data->right_bar_rect, right_bar_height);

    rig_transform_init_identity (data->right_bar_transform);
    rig_transform_translate (data->right_bar_transform,
                             0, data->top_bar_height, 0);

    rig_ui_viewport_set_width (data->tool_vp, data->right_bar_width);
    rig_ui_viewport_set_height (data->tool_vp, right_bar_height);
  }

  rig_rectangle_set_width (data->bottom_bar_rect, data->width - data->right_bar_width);
  rig_rectangle_set_height (data->bottom_bar_rect, data->bottom_bar_height - data->grab_margin);


  rig_transform_init_identity (data->right_bar_transform);
  rig_transform_translate (data->right_bar_transform,
                           data->width - data->right_bar_width + data->grab_margin,
                           data->top_bar_height, 0);

  rig_transform_init_identity (data->main_transform);
  rig_transform_translate (data->main_transform, 0, data->top_bar_height, 0);

  rig_transform_init_identity (data->bottom_bar_transform);
  rig_transform_translate (data->bottom_bar_transform, 0, data->height - data->bottom_bar_height + data->grab_margin, 0);

  //rig_transform_init_identity (data->slider_transform);
  //rig_transform_translate (data->slider_transform, 0, data->bottom_bar_height - 20, 0);

  //rig_transform_init_identity (data->screen_area_transform);
  //rig_transform_translate (data->screen_area_transform, data->main_width / 3, 0, 0);

#if 0
  rig_transform_init_identity (data->pane1_transform);
  rig_transform_translate (data->pane1_transform, data->main_width / 3, 0, 0);

  rig_transform_init_identity (data->pane2_transform);
  rig_transform_translate (data->pane2_transform, (data->main_width / 3) * 2, 0, 0);
#endif

  //rig_slider_set_length (data->slider, data->width);

  rig_arcball_init (&data->arcball,
                    data->main_width / 2,
                    data->main_height / 2,
                    sqrtf (data->main_width * data->main_width + data->main_height * data->main_height) / 2);

  /* picking ray */
  data->picking_ray_color = cogl_pipeline_new (data->ctx->cogl_context);
  cogl_pipeline_set_color4f (data->picking_ray_color, 1.0, 0.0, 0.0, 1.0);
}

static void
data_onscreen_resize (CoglOnscreen *onscreen,
                      int width,
                      int height,
                      void *user_data)
{
  Data *data = user_data;

  data->width = width;
  data->height = height;

  rig_property_dirty (&data->ctx->property_ctx, &data->properties[DATA_PROP_WIDTH]);
  rig_property_dirty (&data->ctx->property_ctx, &data->properties[DATA_PROP_HEIGHT]);

  allocate (data);
}

CoglPipeline *
create_diffuse_specular_material (void)
{
  CoglPipeline *pipeline;
  CoglSnippet *snippet;
  CoglDepthState depth_state;

  pipeline = cogl_pipeline_new (rig_cogl_context);
  cogl_pipeline_set_color4f (pipeline, 0.8f, 0.8f, 0.8f, 1.f);

  /* enable depth testing */
  cogl_depth_state_init (&depth_state);
  cogl_depth_state_set_test_enabled (&depth_state, TRUE);
  cogl_pipeline_set_depth_state (pipeline, &depth_state, NULL);

  /* set up our vertex shader */
  snippet = cogl_snippet_new (COGL_SNIPPET_HOOK_VERTEX,

      /* definitions */
      "uniform mat4 light_shadow_matrix;\n"
      "uniform mat3 normal_matrix;\n"
      "varying vec3 normal_direction, eye_direction;\n"
      "varying vec4 shadow_coords;\n",

      "normal_direction = normalize(normal_matrix * cogl_normal_in);\n"
      "eye_direction    = -vec3(cogl_modelview_matrix * cogl_position_in);\n"

      "shadow_coords = light_shadow_matrix * cogl_modelview_matrix *\n"
      "                cogl_position_in;\n"
);

  cogl_pipeline_add_snippet (pipeline, snippet);
  cogl_object_unref (snippet);

  /* and fragment shader */
  snippet = cogl_snippet_new (COGL_SNIPPET_HOOK_FRAGMENT,
      /* definitions */
      "uniform vec4 light0_ambient, light0_diffuse, light0_specular;\n"
      "uniform vec3 light0_direction_norm;\n"
      "varying vec3 normal_direction, eye_direction;\n",

      /* post */
      NULL);

  cogl_snippet_set_replace (snippet,
      "vec4 final_color = light0_ambient * cogl_color_in;\n"

      " vec3 L = light0_direction_norm;\n"
      " vec3 N = normalize(normal_direction);\n"

      "float lambert = dot(N, L);\n"

      "if (lambert > 0.0)\n"
      "{\n"
      "  final_color += cogl_color_in * light0_diffuse * lambert;\n"

      "  vec3 E = normalize(eye_direction);\n"
      "  vec3 R = reflect (-L, N);\n"
      "  float specular = pow (max(dot(R, E), 0.0),\n"
      "                        2.);\n"
      "  final_color += light0_specular * vec4(.6, .6, .6, 1.0) * specular;\n"
      "}\n"

      "shadow_coords_d = shadow_coords / shadow_coords.w;\n"
      "cogl_texel7 =  cogl_texture_lookup7 (cogl_sampler7, cogl_tex_coord_in[0]);\n"
      "float distance_from_light = cogl_texel7.z + 0.0005;\n"
      "float shadow = 1.0;\n"
      "if (shadow_coords.w > 0.0 && distance_from_light < shadow_coords_d.z)\n"
      "    shadow = 0.5;\n"

      "cogl_color_out = shadow * final_color;\n"
  );

  cogl_pipeline_add_snippet (pipeline, snippet);
  cogl_object_unref (snippet);

  return pipeline;
}

static void
test_init (RigShell *shell, void *user_data)
{
  Data *data = user_data;
  CoglFramebuffer *fb;
  float vector3[3];
  //GError *error = NULL;
  //Item *item;
  //Transition *transition;
  //Path *path;
  int i;
  char *full_path;
  GError *error = NULL;
  CoglPipeline *pipeline;
  RigComponent *component;
  CoglTexture2D *color_buffer;
  CoglPipeline *root_pipeline;
  CoglSnippet *snippet;
  CoglColor color;

  /* A unit test for the list_splice/list_unsplice functions */
#if 0
  _rig_test_list_splice ();
#endif

  data->debug_pick_ray = 1;

  for (i = 0; i < DATA_N_PROPS; i++)
    rig_property_init (&data->properties[i],
                       &data_propert_specs[i],
                       data);

  /* When the items have their own properties then this would
   * be a binding e.g. for the x property of an item. */
  //rig_property_init (&data->path_property,
  //                   &path_property_spec,
  //                   data);

  data->onscreen = cogl_onscreen_new (data->ctx->cogl_context, 880, 660);
  cogl_onscreen_show (data->onscreen);

  /* FIXME: On SDL this isn't taking affect if set before allocating
   * the framebuffer. */
  cogl_onscreen_set_resizable (data->onscreen, TRUE);
  cogl_onscreen_add_resize_handler (data->onscreen, data_onscreen_resize, data);

  fb = COGL_FRAMEBUFFER (data->onscreen);
  data->width = cogl_framebuffer_get_width (fb);
  data->height  = cogl_framebuffer_get_height (fb);


  /*
   * Shadow mapping
   */

  /* Setup the shadow map */
  /* TODO: reallocate if the onscreen framebuffer is resized */
  color_buffer = cogl_texture_2d_new_with_size (rig_cogl_context,
                                                data->width, data->height,
                                                COGL_PIXEL_FORMAT_ANY,
                                                &error);
  if (error)
    g_critical ("could not create texture: %s", error->message);

  data->shadow_color = color_buffer;

  /* XXX: Right now there's no way to disable rendering to the color
   * buffer. */
  data->shadow_fb =
    cogl_offscreen_new_to_texture (COGL_TEXTURE (color_buffer));

  /* retrieve the depth texture */
  cogl_framebuffer_enable_depth_texture (COGL_FRAMEBUFFER (data->shadow_fb),
                                         TRUE);
  /* FIXME: It doesn't seem right that we can query back the texture before
   * the framebuffer has been allocated. */
  data->shadow_map =
    cogl_framebuffer_get_depth_texture (COGL_FRAMEBUFFER (data->shadow_fb));

  if (data->shadow_fb == NULL)
    g_critical ("could not create offscreen buffer");

  /* Hook the shadow sampling */
  root_pipeline = create_diffuse_specular_material ();
  cogl_pipeline_set_layer_texture (root_pipeline, 7, data->shadow_map);

  snippet = cogl_snippet_new (COGL_SNIPPET_HOOK_TEXTURE_LOOKUP,
                              /* declarations */
                              "varying vec4 shadow_coords;\n"
                              "vec4 shadow_coords_d;\n",
                              /* post */
                              "");

  cogl_snippet_set_replace (snippet,
                            "cogl_texel = texture2D(cogl_sampler7, shadow_coords_d.st);\n");

  cogl_pipeline_add_layer_snippet (root_pipeline, 7, snippet);
  cogl_object_unref (snippet);

  data->root_pipeline = root_pipeline;
  data->default_pipeline = cogl_pipeline_new (data->ctx->cogl_context);

  data->diamond_slice_indices =
    cogl_indices_new (data->ctx->cogl_context,
                      COGL_INDICES_TYPE_UNSIGNED_BYTE,
                      _diamond_slice_indices_data,
                      sizeof (_diamond_slice_indices_data) /
                      sizeof (_diamond_slice_indices_data[0]));

  data->circle_texture = rig_create_circle_texture (data->ctx,
                                                    CIRCLE_TEX_RADIUS /* radius */,
                                                    CIRCLE_TEX_PADDING /* padding */);

  data->grid_prim = create_grid (data->ctx,
                                 DEVICE_WIDTH,
                                 DEVICE_HEIGHT,
                                 100,
                                 100);

  data->circle_node_attribute =
    rig_create_circle_fan_p2 (data->ctx, 20, &data->circle_node_n_verts);

  full_path = g_strconcat (RIG_SHARE_DIR, "light-bulb.png", NULL);
  //full_path = g_strconcat (RIG_HANDSET_SHARE_DIR, "nine-slice-test.png", NULL);
  data->light_icon = rig_load_texture (data->ctx, full_path, &error);
  g_free (full_path);
  if (!data->light_icon)
    {
      g_warning ("Failed to load light-bulb texture: %s", error->message);
      g_error_free (error);
    }

  data->timeline_vp = rig_ui_viewport_new (data->ctx, 0, 0, NULL);

  data->screen_area_transform = rig_transform_new (data->ctx, NULL);

  data->device_transform = rig_transform_new (data->ctx, NULL);

  data->camera = rig_camera_new (data->ctx, COGL_FRAMEBUFFER (data->onscreen));
  rig_camera_set_clear (data->camera, FALSE);

  /* XXX: Basically just a hack for now. We should have a
   * RigShellWindow type that internally creates a RigCamera that can
   * be used when handling input events in device coordinates.
   */
  rig_shell_set_window_camera (shell, data->camera);

  data->timeline_camera = rig_camera_new (data->ctx, fb);
  rig_camera_set_clear (data->timeline_camera, FALSE);
  rig_shell_add_input_camera (shell, data->timeline_camera, NULL);
  data->timeline_scale = 1;
  data->timeline_len = 20;

  data->scene = rig_graph_new (data->ctx, NULL);

  data->main_camera_rig0 = rig_entity_new (data->ctx, data->entity_next_id++);
  rig_graphable_add_child (data->scene, data->main_camera_rig0);

  data->main_camera_rig1 = rig_entity_new (data->ctx, data->entity_next_id++);
  rig_graphable_add_child (data->main_camera_rig0, data->main_camera_rig1);

  data->main_camera_rig2 = rig_entity_new (data->ctx, data->entity_next_id++);
  rig_graphable_add_child (data->main_camera_rig1, data->main_camera_rig2);

  data->main_camera_rig3 = rig_entity_new (data->ctx, data->entity_next_id++);
  rig_graphable_add_child (data->main_camera_rig2, data->main_camera_rig3);

  data->origin[0] = DEVICE_WIDTH / 2;
  data->origin[1] = DEVICE_HEIGHT / 2;
  data->origin[2] = 0;

  rig_entity_translate (data->main_camera_rig0,
                        data->origin[0],
                        data->origin[1],
                        data->origin[2]);
                        //DEVICE_WIDTH / 2, (DEVICE_HEIGHT / 2), 0);

  //rig_entity_rotate_z_axis (data->main_camera_rig0, 45);

  rig_entity_translate (data->main_camera_rig2,
                        -DEVICE_WIDTH / 2, -(DEVICE_HEIGHT / 2), 0);

  /* FIXME: currently we also do a z translation due to using
   * cogl_matrix_view_2d_in_perspective, we should stop using that api so we can
   * do our z_2d translation here... */
  data->main_camera_z = 0.f;
  rig_entity_translate (data->main_camera_rig3, 0, 0, data->main_camera_z);

#if 0
  {
    float pos[3] = {0, 0, 0};
    rig_entity_set_position (data->main_camera_rig, pos);
    rig_entity_translate (data->main_camera_rig, 100, 100, 0);
  }
#endif

  data->main_camera = rig_entity_new (data->ctx, data->entity_next_id++);
  //rig_entity_translate (data->main_camera, 100, 100, 0);

#if 0
  data->main_camera_z = 20.f;
  vector3[0] = 0.f;
  vector3[1] = 0.f;
  vector3[2] = data->main_camera_z;
  rig_entity_set_position (data->main_camera, vector3);
#else
  data->main_camera_z = 10.f;
#endif

  data->main_camera_component = rig_camera_new (data->ctx, fb);
  rig_camera_set_clear (data->main_camera_component, FALSE);
  rig_entity_add_component (data->main_camera, data->main_camera_component);
  rig_shell_add_input_camera (shell,
                              data->main_camera_component,
                              data->scene);
                              //data->screen_area_transform);

  data->main_input_region =
    rig_input_region_new_rectangle (0, 0, 0, 0, main_input_region_cb, data);
  rig_camera_add_input_region (data->camera,
                               data->main_input_region);

  rig_graphable_add_child (data->main_camera_rig3, data->main_camera);

  update_camera_position (data);

  /* plane */
  data->plane = rig_entity_new (data->ctx, data->entity_next_id++);
  //data->entities = g_list_prepend (data->entities, data->plane);
  data->pickables = g_list_prepend (data->pickables, data->plane);
  rig_entity_set_cast_shadow (data->plane, FALSE);
  rig_entity_set_y (data->plane, -1.f);

  component = rig_mesh_renderer_new_from_template (data->ctx, "plane");
  rig_entity_add_component (data->plane, component);
  component = rig_material_new_with_pipeline (data->ctx, data->root_pipeline);
  rig_entity_add_component (data->plane, component);

  rig_graphable_add_child (data->scene, data->plane);

  /* 5 cubes */
  pipeline = cogl_pipeline_copy (data->root_pipeline);
  cogl_pipeline_set_color4f (pipeline, 0.6f, 0.6f, 0.6f, 1.0f);
  for (i = 0; i < N_CUBES; i++)
    {
      data->cubes[i] = rig_entity_new (data->ctx, data->entity_next_id++);

      rig_entity_set_scale (data->cubes[i], 50);
      //data->entities = g_list_prepend (data->entities, data->cubes[i]);
      data->pickables = g_list_prepend (data->pickables, data->cubes[i]);

      rig_entity_set_cast_shadow (data->cubes[i], TRUE);
      rig_entity_set_x (data->cubes[i], 50 + i * 125);
      rig_entity_set_y (data->cubes[i], 50 + i * 100);
      rig_entity_set_z (data->cubes[i], 50);
#if 0
      rig_entity_set_y (data->cubes[i], .5);
      rig_entity_set_z (data->cubes[i], 1);
      rig_entity_rotate_y_axis (data->cubes[i], 10);
#endif

      component = rig_mesh_renderer_new_from_template (data->ctx, "cube");
      rig_entity_add_component (data->cubes[i], component);
      component = rig_material_new_with_pipeline (data->ctx, pipeline);
      rig_entity_add_component (data->cubes[i], component);

      rig_graphable_add_child (data->scene, data->cubes[i]);
    }
  cogl_object_unref (pipeline);

  data->light = rig_entity_new (data->ctx, data->entity_next_id++);
  data->entities = g_list_prepend (data->entities, data->light);

  vector3[0] = 0;
  vector3[1] = 0;
  vector3[2] = 200;
  rig_entity_set_position (data->light, vector3);

  rig_entity_rotate_x_axis (data->light, 20);
  rig_entity_rotate_y_axis (data->light, -20);

  component = rig_light_new ();
  cogl_color_init_from_4f (&color, .2f, .2f, .2f, 1.f);
  //cogl_color_init_from_4f (&color, 0, 0, 0, 1.f);
  rig_light_set_ambient (RIG_LIGHT (component), &color);
  cogl_color_init_from_4f (&color, .6f, .6f, .6f, 1.f);
  //cogl_color_init_from_4f (&color, 1, 1, 1, 1.f);
  rig_light_set_diffuse (RIG_LIGHT (component), &color);
  cogl_color_init_from_4f (&color, .4f, .4f, .4f, 1.f);
  //cogl_color_init_from_4f (&color, 1, 1, 1, 1.f);
  rig_light_set_specular (RIG_LIGHT (component), &color);
  rig_light_add_pipeline (RIG_LIGHT (component), root_pipeline);

  rig_entity_add_component (data->light, component);

  component = rig_camera_new (data->ctx, COGL_FRAMEBUFFER (data->shadow_fb));
  data->shadow_map_camera = component;

  rig_camera_set_background_color4f (RIG_CAMERA (component), 0.f, .3f, 0.f, 1.f);
  rig_camera_set_projection_mode (RIG_CAMERA (component),
                                  RIG_PROJECTION_ORTHOGRAPHIC);
  rig_camera_set_orthographic_coordinates (RIG_CAMERA (component),
                                           0, 0, 400, 400);
  rig_camera_set_near_plane (RIG_CAMERA (component), 1.1f);
  rig_camera_set_far_plane (RIG_CAMERA (component), 1000.f);

  rig_entity_add_component (data->light, component);

  rig_graphable_add_child (data->scene, data->light);



  //rig_graphable_add_child (data->main_camera, data->screen_area_transform);
  rig_graphable_add_child (data->screen_area_transform, data->device_transform);

  data->root =
    rig_graph_new (data->ctx,
                   (data->top_bar_transform =
                    rig_transform_new (data->ctx,
                                       (data->top_bar_rect =
                                        rig_rectangle_new4f (data->ctx, 0, 0,
                                                             0.2, 0.2, 0.2, 1)),
                                       NULL)
                   ),
                   (data->left_bar_transform =
                    rig_transform_new (data->ctx,
                                       (data->left_bar_rect =
                                        rig_rectangle_new4f (data->ctx, 0, 0,
                                                             0.2, 0.2, 0.2, 1)),
                                       (data->assets_vp =
                                        rig_ui_viewport_new (data->ctx,
                                                             0, 0,
                                                             NULL)),
                                       NULL)
                   ),
                   (data->right_bar_transform =
                    rig_transform_new (data->ctx,
                                       (data->right_bar_rect =
                                        rig_rectangle_new4f (data->ctx, 0, 0,
                                                             0.2, 0.2, 0.2, 1)),
                                       (data->tool_vp =
                                        rig_ui_viewport_new (data->ctx,
                                                             0, 0,
                                                             NULL)),
                                       NULL)
                   ),
                   (data->main_transform = rig_transform_new (data->ctx, NULL)
                   ),
                   (data->bottom_bar_transform =
                    rig_transform_new (data->ctx,
                                       (data->bottom_bar_rect =
                                        rig_rectangle_new4f (data->ctx, 0, 0,
                                                             0.2, 0.2, 0.2, 1)),
                                       NULL)
                   ),
                   NULL);


  rig_shell_add_input_camera (shell, data->camera, data->root);

#if 0
  data->slider_transform =
    rig_transform_new (data->ctx,
                       data->slider = rig_slider_new (data->ctx,
                                                      RIG_AXIS_X,
                                                      0, 1,
                                                      data->main_width));
  rig_graphable_add_child (data->bottom_bar_transform, data->slider_transform);

  data->slider_progress =
    rig_introspectable_lookup_property (data->slider, "progress");
#endif

  data->timeline_input_region =
    rig_input_region_new_rectangle (0, 0, 0, 0,
                                    timeline_region_input_cb,
                                    data);
  rig_camera_add_input_region (data->timeline_camera,
                               data->timeline_input_region);

  data->timeline = rig_timeline_new (data->ctx, 20.0);
  rig_timeline_set_loop_enabled (data->timeline, TRUE);
  rig_timeline_stop (data->timeline);

  data->timeline_elapsed =
    rig_introspectable_lookup_property (data->timeline, "elapsed");
  data->timeline_progress =
    rig_introspectable_lookup_property (data->timeline, "progress");

#if 0
  rig_property_set_binding (data->slider_progress,
                            update_slider_pos_cb,
                            data,
                            data->timeline_elapsed,
                            NULL);

  rig_property_set_binding (data->timeline_progress,
                            update_timeline_progress_cb,
                            data,
                            data->slider_progress,
                            NULL);
#endif

  //rig_graphable_add_child (data->camera, RIG_OBJECT (data->root));

#if 0
  item = item_new (data, "contact0.png", 0);
  rig_graphable_add_child (data->device_transform, item->transform);
  data->items = g_list_prepend (data->items, item);

  item = item_new (data, "contact1.png", 1);
  data->items = g_list_prepend (data->items, item);
  rig_graphable_add_child (data->device_transform, item->transform);

  data->selected_item = item;

  transition = transition_new (data, 0);
  data->transitions = g_list_prepend (data->transitions, transition);
  data->selected_transition = transition;

  path = path_new_for_property (data->ctx,
                                &transition->props[TRANSITION_PROP_PROGRESS],
                                &item->props[ITEM_PROP_X]);
  path_insert_float (path, -1, 0);
  path_insert_float (path, 0, 0);
  path_insert_float (path, 5, 100);
  path_insert_float (path, 10, 200);
  path_insert_float (path, 20, 200);
  path_insert_float (path, 21, 200);
  transition_add_path (transition, path);

  path = path_new_for_property (data->ctx,
                                &transition->props[TRANSITION_PROP_PROGRESS],
                                &item->props[ITEM_PROP_Y]);
  path_insert_float (path, -1, 0);
  path_insert_float (path, 0, 0);
  path_insert_float (path, 5, 100);
  path_insert_float (path, 10, 0);
  path_insert_float (path, 20, 0);
  path_insert_float (path, 21, 0);
  transition_add_path (transition, path);

  path = path_new_for_property (data->ctx,
                                &transition->props[TRANSITION_PROP_PROGRESS],
                                &item->props[ITEM_PROP_Z]);
  path_insert_float (path, -1, 0);
  path_insert_float (path, 0, 0);
  path_insert_float (path, 5, 0);
  path_insert_float (path, 10, 0);
  path_insert_float (path, 20, 0);
  path_insert_float (path, 21, 0);
  transition_add_path (transition, path);

  path = path_new_for_property (data->ctx,
                                &transition->props[TRANSITION_PROP_PROGRESS],
                                &item->props[ITEM_PROP_ROTATION]);
  path_insert_quaternion (path, -1, 0, 0, 1, 0);
  path_insert_quaternion (path, 0, 0, 0, 1, 0);
  path_insert_quaternion (path, 10, 90, 0, 1, 0);
  path_insert_quaternion (path, 20, 0, 0, 1, 0);
  path_insert_quaternion (path, 21, 0, 0, 1, 0);
  transition_add_path (transition, path);
#endif

    /* tool */
  data->tool = rig_tool_new (data->shell);
  rig_tool_set_camera (data->tool, data->main_camera);

  allocate (data);

#ifndef __ANDROID__
  if (_rig_handset_remaining_args &&
      _rig_handset_remaining_args[0])
    {
      char *ui;
      struct stat st;

      stat (_rig_handset_remaining_args[0], &st);
      if (!S_ISDIR (st.st_mode))
        {
          g_error ("Could not find project directory %s",
                   _rig_handset_remaining_args[0]);
        }

      _rig_project_dir = _rig_handset_remaining_args[0];

      ui = g_strconcat (_rig_handset_remaining_args[0], "/ui.xml");

      load (data, ui);
      g_free (ui);
    }
#endif
}

static void
test_fini (RigShell *shell, void *user_data)
{
  Data *data = user_data;
  int i;

  rig_ref_countable_unref (data->camera);
  rig_ref_countable_unref (data->root);

  for (i = 0; i < DATA_N_PROPS; i++)
    rig_property_destroy (&data->properties[i]);

  rig_ref_countable_unref (data->timeline_vp);

  cogl_object_unref (data->circle_texture);

  cogl_object_unref (data->grid_prim);
  cogl_object_unref (data->circle_node_attribute);

  cogl_object_unref (data->light_icon);

  //cogl_object_unref (data->rounded_tex);
}

static RigInputEventStatus
test_input_handler (RigInputEvent *event, void *user_data)
{
  Data *data = user_data;

  if (rig_input_event_get_type (event) == RIG_INPUT_EVENT_TYPE_MOTION)
    {
      /* Anything that can claim the keyboard focus will do so during
       * motion events so we clear it before running other input
       * callbacks */
      data->key_focus_callback = NULL;
    }

  switch (rig_input_event_get_type (event))
    {
    case RIG_INPUT_EVENT_TYPE_MOTION:
#if 0
      switch (rig_motion_event_get_action (event))
        {
        case RIG_MOTION_EVENT_ACTION_DOWN:
          //g_print ("Press Down\n");
          return RIG_INPUT_EVENT_STATUS_HANDLED;
        case RIG_MOTION_EVENT_ACTION_UP:
          //g_print ("Release\n");
          return RIG_INPUT_EVENT_STATUS_HANDLED;
        case RIG_MOTION_EVENT_ACTION_MOVE:
          //g_print ("Move\n");
          return RIG_INPUT_EVENT_STATUS_HANDLED;
        }
#endif
      break;

    case RIG_INPUT_EVENT_TYPE_KEY:
      {
        if (data->key_focus_callback)
          data->key_focus_callback (event, data);
      }
      break;
    }

  return RIG_INPUT_EVENT_STATUS_UNHANDLED;
}

typedef struct _SaveState
{
  Data *data;
  FILE *file;
  int indent;
  RigEntity *current_entity;
} SaveState;

static void
save_component_cb (RigComponent *component,
                   void *user_data)
{
  const RigType *type = rig_object_get_type (component);
  SaveState *state = user_data;

  state->indent += INDENT_LEVEL;

  if (type == &rig_light_type)
    {
      RigLight *light = RIG_LIGHT (component);
      const CoglColor *ambient = rig_light_get_ambient (light);
      const CoglColor *diffuse = rig_light_get_diffuse (light);
      const CoglColor *specular = rig_light_get_specular (light);

      fprintf (state->file,
               "%*s<light "
               "ambient=\"#%02x%02x%02x%02x\" "
               "diffuse=\"#%02x%02x%02x%02x\" "
               "specular\"#%02x%02x%02x%02x\"/>\n",
               state->indent, "",
               cogl_color_get_red_byte (ambient),
               cogl_color_get_green_byte (ambient),
               cogl_color_get_blue_byte (ambient),
               cogl_color_get_alpha_byte (ambient),
               cogl_color_get_red_byte (diffuse),
               cogl_color_get_green_byte (diffuse),
               cogl_color_get_blue_byte (diffuse),
               cogl_color_get_alpha_byte (diffuse),
               cogl_color_get_red_byte (specular),
               cogl_color_get_green_byte (specular),
               cogl_color_get_blue_byte (specular),
               cogl_color_get_alpha_byte (specular));
    }
  else if (type == &rig_material_type)
    {
      /* FIXME: The RigMaterial should maintain a reference to the original
       * asset so we can get its id here. */
      g_warning ("FIXME: Unable to save materials correctly\n");
      fprintf (state->file, "%*s<material><texture asset=\"0\"/></material>\n",
               state->indent, "");
    }
  else if (type == &rig_diamond_type)
    {
      fprintf (state->file, "%*s<diamond size=\"%f\"/>\n",
               state->indent, "",
               rig_diamond_get_size (RIG_DIAMOND (component)));
    }

  state->indent -= INDENT_LEVEL;
}

static RigTraverseVisitFlags
_rig_entitygraph_pre_save_cb (RigObject *object,
                              int depth,
                              void *user_data)
{
  SaveState *state = user_data;
  RigType *type = rig_object_get_type (object);
  RigObject *parent = rig_graphable_get_parent (object);
  RigEntity *entity;
  CoglQuaternion *q;
  float angle;
  float axis[3];
  GList *l2;

  if (type != &rig_entity_type)
    {
      g_warning ("Can't save non-entity graphables\n");
      return RIG_TRAVERSE_VISIT_CONTINUE;
    }

  entity = object;

  state->indent += INDENT_LEVEL;
  fprintf (state->file, "%*s<entity id=\"%d\"",
           state->indent, "",
           rig_entity_get_id (entity));

  if (parent && rig_object_get_type (parent) == &rig_entity_type)
    {
      RigEntity *parent_entity = parent;
      fprintf (state->file, " parent=\"%d\"", rig_entity_get_id (parent));
    }

  q = rig_entity_get_rotation (entity);

  angle = cogl_quaternion_get_rotation_angle (q);
  cogl_quaternion_get_rotation_axis (q, axis);

  fprintf (state->file,
           " x=\"%f\" y=\"%f\" z=\"%f\" "
           "rotation=\"[%f (%f, %f, %f)]]\">\n",
           rig_entity_get_x (entity),
           rig_entity_get_y (entity),
           rig_entity_get_z (entity),
           angle, axis[0], axis[1], axis[2]);

  state->current_entity = entity;
  rig_entity_foreach_component (entity,
                                save_component_cb,
                                state);

  fprintf (state->file, "%*s</entity>\n", state->indent, "");
  state->indent -= INDENT_LEVEL;

  return RIG_TRAVERSE_VISIT_CONTINUE;
}

static void
save (Data *data)
{
  const char *project_name = "Flibble";
  struct stat sb;
  char *path = g_strdup_printf ("%s/ui.xml", project_name);
  FILE *file;
  SaveState state;
  GList *l;

  if (stat (project_name, &sb) == -1)
    mkdir (project_name, 0777);

  file = fopen (path, "w");
  if (!file)
    {
      g_warning ("Failed to open %s file for saving\n", path);
      return;
    }

  state.data = data;
  state.file = file;
  state.indent = 0;

  fprintf (file, "<ui>\n");

  /* Assets */

  for (l = data->assets; l; l = l->next)
    {
      Asset *asset = l->data;

      if (asset->type != ASSET_TYPE_TEXTURE)
        continue;

      state.indent += INDENT_LEVEL;
      fprintf (file, "%*s<asset id=\"%d\" type=\"texture\" path=\"%s\" />\n",
               state.indent, "",
               asset->id,
               asset->path);
    }

  for (l = data->items; l; l = l->next)
    {
      Item *item = l->data;
      //GList *l2;
      //int i;

      /* TODO: Items should reference assets */

      state.indent += INDENT_LEVEL;
      fprintf (file, "%*s<item id=\"%d\">\n", state.indent, "", item->id);

      state.indent += INDENT_LEVEL;
      fprintf (file, "%*s<texture asset=\"%d\"/>\n", state.indent, "", item->texture_asset->id);
      state.indent -= INDENT_LEVEL;

      fprintf (file, "%*s</item>\n", state.indent, "");
      state.indent -= INDENT_LEVEL;
    }

  rig_graphable_traverse (data->scene,
                          RIG_TRAVERSE_DEPTH_FIRST,
                          _rig_entitygraph_pre_save_cb,
                          NULL,
                          &state);

  for (l = data->transitions; l; l = l->next)
    {
      Transition *transition = l->data;
      GList *l2;
      //int i;

      state.indent += INDENT_LEVEL;
      fprintf (file, "%*s<transition id=\"%d\">\n", state.indent, "", transition->id);

      for (l2 = transition->paths; l2; l2 = l2->next)
        {
          Path *path = l2->data;
          GList *l3;
          Item *item;

          if (path == NULL)
            continue;

          item = path->prop->object;

          state.indent += INDENT_LEVEL;
          fprintf (file, "%*s<path item=\"%d\" property=\"%s\">\n", state.indent, "",
                   item->id,
                   path->prop->spec->name);

          state.indent += INDENT_LEVEL;
          for (l3 = path->nodes.head; l3; l3 = l3->next)
            {
              switch (path->prop->spec->type)
                {
                case RIG_PROPERTY_TYPE_FLOAT:
                  {
                    NodeFloat *node = l3->data;
                    fprintf (file, "%*s<node t=\"%f\" value=\"%f\" />\n", state.indent, "", node->t, node->value);
                    break;
                  }
                case RIG_PROPERTY_TYPE_QUATERNION:
                  {
                    NodeQuaternion *node = l3->data;
                    CoglQuaternion *q = &node->value;
                    float angle;
                    float axis[3];

                    angle = cogl_quaternion_get_rotation_angle (q);
                    cogl_quaternion_get_rotation_axis (q, axis);

                    fprintf (file, "%*s<node t=\"%f\" value=\"[%f (%f, %f, %f)]\" />\n", state.indent, "",
                             node->t, angle, axis[0], axis[1], axis[2]);
                    break;
                  }
                default:
                  g_warn_if_reached ();
                }
            }
          state.indent -= INDENT_LEVEL;

          fprintf (file, "%*s</path>\n", state.indent, "");
          state.indent -= INDENT_LEVEL;
        }

      fprintf (file, "%*s</transition>\n", state.indent, "");
      fprintf (file, "</ui>\n");
    }

  fclose (file);

  g_print ("File Saved\n");
}

static RigInputEventStatus
asset_input_cb (RigInputRegion *region,
                RigInputEvent *event,
                void *user_data)
{
  Asset *asset = user_data;
  Data *data = asset->data;

  g_print ("Asset input\n");

  if (rig_input_event_get_type (event) == RIG_INPUT_EVENT_TYPE_MOTION)
    {
      if (rig_motion_event_get_action (event) == RIG_MOTION_EVENT_ACTION_DOWN)
        {
          Item *item = item_new (data,
                                 asset,
                                 data->items_next_id++);

          RigEntity *entity = rig_entity_new (data->ctx,
                                              data->entity_next_id++);
          RigComponent *component =
            rig_material_new_with_texture (data->ctx,
                                           asset->texture);
          rig_entity_add_component (entity, component);
          component = rig_diamond_new (data->ctx,
                                       400,
                                       cogl_texture_get_width (asset->texture),
                                       cogl_texture_get_height (asset->texture));
          rig_entity_add_component (entity, component);

          //data->entities = g_list_prepend (data->entities, entity);
          data->selected_entity = entity;
          //rig_graphable_add_child (data->device_transform, entity);
          rig_graphable_add_child (data->scene, entity);

          data->items = g_list_prepend (data->items, item);
          data->selected_item = item;
          rig_graphable_add_child (data->device_transform, item);
          rig_shell_queue_redraw (data->ctx->shell);
          return RIG_INPUT_EVENT_STATUS_HANDLED;
        }
    }

  return RIG_INPUT_EVENT_STATUS_UNHANDLED;
}

static RigInputEventStatus
add_light_cb (RigInputRegion *region,
              RigInputEvent *event,
              void *user_data)
{
  if (rig_input_event_get_type (event) == RIG_INPUT_EVENT_TYPE_MOTION)
    {
      if (rig_motion_event_get_action (event) == RIG_MOTION_EVENT_ACTION_DOWN)
        {
          g_print ("Add light!\n");
          return RIG_INPUT_EVENT_STATUS_HANDLED;
        }
    }

  return RIG_INPUT_EVENT_STATUS_UNHANDLED;
}

static void
add_asset_icon (Data *data,
                CoglTexture *texture,
                float y_pos,
                RigInputRegionCallback callback,
                void *user_data)
{
  RigNineSlice *nine_slice;
  RigInputRegion *region;
  RigTransform *transform =
    rig_transform_new (data->ctx,
                       (nine_slice = rig_nine_slice_new (data->ctx,
                                                         texture,
                                                         0, 0, 0, 0,
                                                         100, 100)),
                       (region =
                        rig_input_region_new_rectangle (0, 0, 100, 100,
                                                        callback,
                                                        user_data)),
                       NULL);
  rig_graphable_add_child (data->assets_list, transform);

  rig_transform_translate (transform, 10, y_pos, 0);

  //rig_input_region_set_graphable (region, nine_slice);

  rig_ref_countable_unref (transform);
  rig_ref_countable_unref (nine_slice);
  rig_ref_countable_unref (region);

}

static void
update_asset_list (Data *data)
{
  GList *l;
  int i = 0;
  RigObject *doc_node;

  if (data->assets_list)
    rig_graphable_remove_child (data->assets_list);

  data->assets_list = rig_graph_new (data->ctx, NULL);

  doc_node = rig_ui_viewport_get_doc_node (data->assets_vp);
  rig_graphable_add_child (doc_node, data->assets_list);
  rig_ref_countable_unref (data->assets_list);

  for (l = data->assets, i= 0; l; l = l->next, i++)
    {
      Asset *asset = l->data;
      add_asset_icon (data, asset->texture, 10 + 110 * i, asset_input_cb, asset);
    }

  add_asset_icon (data, data->light_icon, 10 + 110 * i++, add_light_cb, NULL);
}

enum {
  LOADER_STATE_NONE,
  LOADER_STATE_LOADING_ITEM,
  LOADER_STATE_LOADING_ENTITY,
  LOADER_STATE_LOADING_MATERIAL_COMPONENT,
  LOADER_STATE_LOADING_MESH_COMPONENT,
  LOADER_STATE_LOADING_DIAMOND_COMPONENT,
  LOADER_STATE_LOADING_LIGHT_COMPONENT,
  LOADER_STATE_LOADING_CAMERA_COMPONENT,
  LOADER_STATE_LOADING_TRANSITION,
  LOADER_STATE_LOADING_PATH
};

typedef struct _Loader
{
  Data *data;
  GQueue state;
  uint32_t id;
  CoglBool texture_specified;
  uint32_t texture_asset_id;

  GList *assets;
  GList *items;
  GList *entities;
  GList *lights;
  GList *transitions;


  float diamond_size;
  RigEntity *current_entity;
  CoglBool is_light;

  Transition *current_transition;
  Path *current_path;

  GHashTable *id_map;
} Loader;

static void
loader_push_state (Loader *loader, int state)
{
  g_queue_push_tail (&loader->state, GINT_TO_POINTER (state));
}

static int
loader_get_state (Loader *loader)
{
  void *state = g_queue_peek_tail (&loader->state);
  return GPOINTER_TO_INT (state);
}

static void
loader_pop_state (Loader *loader)
{
  g_queue_pop_tail (&loader->state);
}

static Item *
loader_find_item (Loader *loader, uint32_t id)
{
  GList *l;

  for (l = loader->items; l; l = l->next)
    {
      Item *item = l->data;

      if (item->id == id)
        return item;
    }

  return NULL;
}

static Asset *
loader_find_asset (Loader *loader, uint32_t id)
{
  GList *l;

  for (l = loader->assets; l; l = l->next)
    {
      Asset *asset = l->data;

      if (asset->id == id)
        return asset;
    }

  return NULL;
}


static void
parse_start_element (GMarkupParseContext *context,
                     const char *element_name,
                     const char **attribute_names,
                     const char **attribute_values,
                     void *user_data,
                     GError **error)
{
  Loader *loader = user_data;
  Data *data = loader->data;
  int state = loader_get_state (loader);

  if (state == LOADER_STATE_NONE &&
      strcmp (element_name, "asset") == 0)
    {
      const char *id_str;
      const char *type;
      const char *path;
      uint32_t id;

      if (!g_markup_collect_attributes (element_name,
                                        attribute_names,
                                        attribute_values,
                                        error,
                                        G_MARKUP_COLLECT_STRING,
                                        "id",
                                        &id_str,
                                        G_MARKUP_COLLECT_STRING,
                                        "type",
                                        &type,
                                        G_MARKUP_COLLECT_STRING,
                                        "path",
                                        &path,
                                        G_MARKUP_COLLECT_INVALID))
        {
          return;
        }

      id = g_ascii_strtoull (id_str, NULL, 10);

      if (strcmp (type, "texture") == 0)
        {
          Asset *asset = asset_new_texture (data, id, path);
          loader->assets = g_list_prepend (loader->assets, asset);
        }
      else
        g_warning ("Ignoring unknown asset type: %s\n", type);
    }
  else if (state == LOADER_STATE_NONE &&
           strcmp (element_name, "item") == 0)
    {
      const char *id_str;

      loader_push_state (loader, LOADER_STATE_LOADING_ITEM);

      if (!g_markup_collect_attributes (element_name,
                                        attribute_names,
                                        attribute_values,
                                        error,
                                        G_MARKUP_COLLECT_STRING,
                                        "id",
                                        &id_str,
                                        G_MARKUP_COLLECT_INVALID))
        {
          return;
        }

      loader->id = g_ascii_strtoull (id_str, NULL, 10);
    }
  else if (state == LOADER_STATE_NONE &&
           strcmp (element_name, "entity") == 0)
    {
      const char *id_str;
      unsigned int id;
      RigEntity *entity;
      const char *parent_id_str;
      const char *x_str;
      const char *y_str;
      const char *z_str;
      const char *rotation_str;
      const char *scale_str;

      if (!g_markup_collect_attributes (element_name,
                                        attribute_names,
                                        attribute_values,
                                        error,
                                        G_MARKUP_COLLECT_STRING,
                                        "id",
                                        &id_str,
                                        G_MARKUP_COLLECT_STRING|G_MARKUP_COLLECT_OPTIONAL,
                                        "parent",
                                        &parent_id_str,
                                        G_MARKUP_COLLECT_STRING|G_MARKUP_COLLECT_OPTIONAL,
                                        "x",
                                        &x_str,
                                        G_MARKUP_COLLECT_STRING|G_MARKUP_COLLECT_OPTIONAL,
                                        "y",
                                        &y_str,
                                        G_MARKUP_COLLECT_STRING|G_MARKUP_COLLECT_OPTIONAL,
                                        "z",
                                        &z_str,
                                        G_MARKUP_COLLECT_STRING|G_MARKUP_COLLECT_OPTIONAL,
                                        "rotation",
                                        &rotation_str,
                                        G_MARKUP_COLLECT_STRING|G_MARKUP_COLLECT_OPTIONAL,
                                        "scale",
                                        &scale_str,
                                        G_MARKUP_COLLECT_INVALID))
      {
        return;
      }

      id = g_ascii_strtoull (id_str, NULL, 10);
      if (g_hash_table_lookup (loader->id_map, GUINT_TO_POINTER (id)))
        {
          g_set_error (error,
                       G_MARKUP_ERROR,
                       G_MARKUP_ERROR_INVALID_CONTENT,
                       "Duplicate entity id %d", id);
          return;
        }

      entity = rig_entity_new (loader->data->ctx, loader->data->entity_next_id++);

      if (parent_id_str)
        {
          unsigned int parent_id = g_ascii_strtoull (parent_id_str, NULL, 10);
          RigEntity *parent = g_hash_table_lookup (loader->id_map,
                                                   GUINT_TO_POINTER (parent_id));
          if (!parent)
            {
              g_set_error (error,
                           G_MARKUP_ERROR,
                           G_MARKUP_ERROR_INVALID_CONTENT,
                           "Invalid parent id referenced in entity element");
              rig_ref_countable_unref (entity);
              return;
            }

          rig_graphable_add_child (parent, entity);
        }

      if (x_str)
        {
          double x = g_ascii_strtod (x_str, NULL);
          rig_entity_set_x (entity, x);
        }
      if (y_str)
        {
          double y = g_ascii_strtod (y_str, NULL);
          rig_entity_set_y (entity, y);
        }
      if (z_str)
        {
          double z = g_ascii_strtod (z_str, NULL);
          rig_entity_set_z (entity, z);
        }
      if (rotation_str)
        {
          float angle;
          float axis[3];
          CoglQuaternion q;

          if (sscanf (rotation_str, "[%f (%f, %f, %f)]",
                      &angle, &axis[0], &axis[1], &axis[2]) != 4)
            {
              g_set_error (error,
                           G_MARKUP_ERROR,
                           G_MARKUP_ERROR_INVALID_CONTENT,
                           "Invalid entity rotation");
              return;
            }
          cogl_quaternion_init_from_angle_vector (&q, angle, axis);
          rig_entity_set_rotation (entity, &q);
        }
      if (scale_str)
        {
          double scale = g_ascii_strtod (scale_str, NULL);
          rig_entity_set_scale (entity, scale);
        }

      loader->current_entity = entity;
      loader->is_light = FALSE;
      g_hash_table_insert (loader->id_map, GUINT_TO_POINTER (id), entity);

      loader_push_state (loader, LOADER_STATE_LOADING_ENTITY);
    }
  else if (state == LOADER_STATE_LOADING_ENTITY &&
           strcmp (element_name, "material") == 0)
    {
      loader->texture_specified = FALSE;
      loader_push_state (loader, LOADER_STATE_LOADING_MATERIAL_COMPONENT);
#if 0
      const char *texture_id;

      g_markup_collect_attributes (element_name,
                                   attribute_names,
                                   attribute_values,
                                   error,
                                   G_MARKUP_COLLECT_STRING|G_MARKUP_COLLECT_OPTIONAL,
                                   "texture",
                                   &texture_id,
                                   G_MARKUP_COLLECT_INVALID);
#endif
    }
  else if (state == LOADER_STATE_LOADING_ENTITY &&
           strcmp (element_name, "light") == 0)
    {
      const char *ambient_str;
      const char *diffuse_str;
      const char *specular_str;
      CoglColor ambient;
      CoglColor diffuse;
      CoglColor specular;
      RigLight *light;

      if (!g_markup_collect_attributes (element_name,
                                        attribute_names,
                                        attribute_values,
                                        error,
                                        G_MARKUP_COLLECT_STRING,
                                        "ambient",
                                        &ambient_str,
                                        G_MARKUP_COLLECT_STRING,
                                        "diffuse",
                                        &diffuse_str,
                                        G_MARKUP_COLLECT_STRING,
                                        "specular",
                                        &specular_str,
                                        G_MARKUP_COLLECT_INVALID))
        {
          return;
        }

      rig_util_parse_color (loader->data->ctx, ambient_str, &ambient);
      rig_util_parse_color (loader->data->ctx, diffuse_str, &diffuse);
      rig_util_parse_color (loader->data->ctx, specular_str, &specular);

      light = rig_light_new ();
      rig_light_set_ambient (light, &ambient);
      rig_light_set_diffuse (light, &diffuse);
      rig_light_set_specular (light, &specular);

      rig_entity_add_component (loader->current_entity, light);
      loader->is_light = TRUE;

      //loader_push_state (loader, LOADER_STATE_LOADING_LIGHT_COMPONENT);
    }
  else if (state == LOADER_STATE_LOADING_ENTITY &&
           strcmp (element_name, "diamond") == 0)
    {
      const char *size_str;

      if (!g_markup_collect_attributes (element_name,
                                        attribute_names,
                                        attribute_values,
                                        error,
                                        G_MARKUP_COLLECT_STRING,
                                        "size",
                                        &size_str,
                                        G_MARKUP_COLLECT_INVALID))
        return;

      loader->diamond_size = g_ascii_strtod (size_str, NULL);

      loader_push_state (loader, LOADER_STATE_LOADING_DIAMOND_COMPONENT);
    }
  else if ((state == LOADER_STATE_LOADING_ITEM ||
            state == LOADER_STATE_LOADING_MATERIAL_COMPONENT) &&
           strcmp (element_name, "texture") == 0)
    {
      const char *id_str;

      if (!g_markup_collect_attributes (element_name,
                                        attribute_names,
                                        attribute_values,
                                        error,
                                        G_MARKUP_COLLECT_STRING,
                                        "asset",
                                        &id_str,
                                        G_MARKUP_COLLECT_INVALID))
        return;

      loader->texture_specified = TRUE;
      loader->texture_asset_id = g_ascii_strtoull (id_str, NULL, 10);
    }
  else if (state == LOADER_STATE_NONE &&
           strcmp (element_name, "transition") == 0)
    {
      const char *id_str;
      uint32_t id;

      loader_push_state (loader, LOADER_STATE_LOADING_TRANSITION);

      if (!g_markup_collect_attributes (element_name,
                                        attribute_names,
                                        attribute_values,
                                        error,
                                        G_MARKUP_COLLECT_STRING,
                                        "id",
                                        &id_str,
                                        G_MARKUP_COLLECT_INVALID))
        return;

      id = g_ascii_strtoull (id_str, NULL, 10);

      loader->current_transition = transition_new (loader->data, id);
      loader->transitions = g_list_prepend (loader->transitions, loader->current_transition);
    }
  else if (state == LOADER_STATE_LOADING_TRANSITION &&
           strcmp (element_name, "path") == 0)
    {
      const char *item_id_str;
      uint32_t item_id;
      Item *item;
      const char *property_name;
      RigProperty *prop;

      if (!g_markup_collect_attributes (element_name,
                                        attribute_names,
                                        attribute_values,
                                        error,
                                        G_MARKUP_COLLECT_STRING,
                                        "item",
                                        &item_id_str,
                                        G_MARKUP_COLLECT_STRING,
                                        "property",
                                        &property_name,
                                        G_MARKUP_COLLECT_INVALID))
        return;

      item_id = g_ascii_strtoull (item_id_str, NULL, 10);

      item = loader_find_item (loader, item_id);
      if (!item)
        {
          g_set_error (error,
                       G_MARKUP_ERROR,
                       G_MARKUP_ERROR_INVALID_CONTENT,
                       "Invalid Item id referenced in path element");
          return;
        }

      prop = rig_introspectable_lookup_property (item, property_name);
      if (!prop)
        {
          g_set_error (error,
                       G_MARKUP_ERROR,
                       G_MARKUP_ERROR_INVALID_CONTENT,
                       "Invalid Item property referenced in path element");
          return;
        }

      loader->current_path =
        path_new_for_property (data->ctx,
                               &loader->current_transition->props[TRANSITION_PROP_PROGRESS],
                               prop);

      loader_push_state (loader, LOADER_STATE_LOADING_PATH);
    }
  else if (state == LOADER_STATE_LOADING_PATH &&
           strcmp (element_name, "node") == 0)
    {
      const char *t_str;
      float t;
      const char *value_str;

      if (!g_markup_collect_attributes (element_name,
                                        attribute_names,
                                        attribute_values,
                                        error,
                                        G_MARKUP_COLLECT_STRING,
                                        "t",
                                        &t_str,
                                        G_MARKUP_COLLECT_STRING,
                                        "value",
                                        &value_str,
                                        G_MARKUP_COLLECT_INVALID))
        return;

      t = g_ascii_strtod (t_str, NULL);

      switch (loader->current_path->prop->spec->type)
        {
        case RIG_PROPERTY_TYPE_FLOAT:
          {
            float value = g_ascii_strtod (value_str, NULL);
            path_insert_float (loader->current_path, t, value);
            break;
          }
        case RIG_PROPERTY_TYPE_QUATERNION:
          {
            float angle, x, y, z;

            if (sscanf (value_str, "[%f (%f, %f, %f)]", &angle, &x, &y, &z) != 4)
              {
                g_set_error (error,
                             G_MARKUP_ERROR,
                             G_MARKUP_ERROR_INVALID_CONTENT,
                             "Invalid rotation value");
                return;
              }

            path_insert_quaternion (loader->current_path, t, angle, x, y, z);
            break;
          }
        }
    }
}

static void
parse_end_element (GMarkupParseContext *context,
                   const char *element_name,
                   void *user_data,
                   GError **error)
{
  Loader *loader = user_data;
  int state = loader_get_state (loader);

  if (state == LOADER_STATE_LOADING_ITEM &&
      strcmp (element_name, "item") == 0)
    {
      Item *item;
      //Data *data = loader->data;
      Asset *texture_asset;
      RigEntity *entity;
      RigComponent *component;

      texture_asset = loader_find_asset (loader, loader->texture_asset_id);
      if (!texture_asset)
        {
          g_set_error (error,
                       G_MARKUP_ERROR,
                       G_MARKUP_ERROR_INVALID_CONTENT,
                       "Invalid asset id");
          return;
        }

      item = item_new (loader->data,
                       texture_asset,
                       loader->id);
      entity = rig_entity_new (loader->data->ctx,
                               loader->data->entity_next_id++);
      component = rig_material_new_with_texture (loader->data->ctx,
                                                 texture_asset->texture);
      rig_entity_add_component (entity, component);
      component = rig_diamond_new (loader->data->ctx,
                                   400,
                                   cogl_texture_get_width (texture_asset->texture),
                                   cogl_texture_get_height (texture_asset->texture));
      rig_entity_add_component (entity, component);

      g_print ("loaded item = %p\n", item);

      loader->items = g_list_prepend (loader->items, item);
      loader->entities = g_list_prepend (loader->entities, entity);

      loader_pop_state (loader);
    }
  else if (state == LOADER_STATE_LOADING_ENTITY &&
           strcmp (element_name, "entity") == 0)
    {
      if (loader->is_light)
        loader->lights = g_list_prepend (loader->entities, loader->current_entity);
      else
        loader->entities = g_list_prepend (loader->entities, loader->current_entity);

      loader_pop_state (loader);
    }
  else if (state == LOADER_STATE_LOADING_DIAMOND_COMPONENT &&
           strcmp (element_name, "diamond") == 0)
    {
      RigMaterial *material =
        rig_entity_get_component (loader->current_entity,
                                  RIG_COMPONENT_TYPE_MATERIAL);
      CoglTexture *texture = NULL;
      RigDiamond *diamond;

      /* We need to know the size of the texture before we can create
       * a diamond component */
      if (material)
        {
          CoglPipeline *pipeline = rig_material_get_pipeline (material);
          texture = cogl_pipeline_get_layer_texture (pipeline, 0);
        }

      if (!texture)
        {
          g_set_error (error,
                       G_MARKUP_ERROR,
                       G_MARKUP_ERROR_INVALID_CONTENT,
                       "Can't add diamond component without a texture");
          return;
        }

      diamond = rig_diamond_new (loader->data->ctx,
                                 loader->diamond_size,
                                 cogl_texture_get_width (texture),
                                 cogl_texture_get_height (texture));
      rig_entity_add_component (loader->current_entity,
                                diamond);

      loader_pop_state (loader);
    }
  else if (state == LOADER_STATE_LOADING_MATERIAL_COMPONENT &&
           strcmp (element_name, "material") == 0)
    {
      CoglPipeline *pipeline = cogl_pipeline_new (loader->data->ctx->cogl_context);
      RigMaterial *material;
      Asset *texture_asset;

      if (loader->texture_specified)
        {
          texture_asset = loader_find_asset (loader, loader->texture_asset_id);
          if (!texture_asset)
            {
              g_set_error (error,
                           G_MARKUP_ERROR,
                           G_MARKUP_ERROR_INVALID_CONTENT,
                           "Invalid asset id");
              return;
            }

          cogl_pipeline_set_layer_texture (pipeline, 0, texture_asset->texture);
        }

      material = rig_material_new_with_pipeline (loader->data->ctx, pipeline);
      rig_entity_add_component (loader->current_entity, material);

      loader_pop_state (loader);
    }
  else if (state == LOADER_STATE_LOADING_TRANSITION &&
           strcmp (element_name, "transition") == 0)
    {
      loader_pop_state (loader);
    }
  else if (state == LOADER_STATE_LOADING_PATH &&
           strcmp (element_name, "path") == 0)
    {
      transition_add_path (loader->current_transition, loader->current_path);
      loader_pop_state (loader);
    }
}

static void
parse_error (GMarkupParseContext *context,
             GError *error,
             void *user_data)
{

}

static void
free_ux (Data *data)
{
  GList *l;

  for (l = data->transitions; l; l = l->next)
    transition_free (l->data);
  g_list_free (data->transitions);
  data->transitions = NULL;

  for (l = data->items; l; l = l->next)
    item_free (l->data);
  g_list_free (data->items);
  data->items = NULL;

  for (l = data->assets; l; l = l->next)
    _asset_free (l->data);
  g_list_free (data->assets);
  data->assets = NULL;
}

static void
load (Data *data, const char *file)
{
  GMarkupParser parser = {
    .start_element = parse_start_element,
    .end_element = parse_end_element,
    .error = parse_error
  };
  Loader loader;
  char *contents;
  gsize len;
  GError *error = NULL;
  GMarkupParseContext *context;
  GList *l;

  memset (&loader, 0, sizeof (loader));
  loader.data = data;
  g_queue_init (&loader.state);
  loader_push_state (&loader, LOADER_STATE_NONE);
  g_queue_push_tail (&loader.state, LOADER_STATE_NONE);

  loader.id_map = g_hash_table_new (g_direct_hash, g_direct_equal);

  if (!g_file_get_contents (file,
                            &contents,
                            &len,
                            &error))
    {
      g_warning ("Failed to load ui description: %s", error->message);
      return;
    }

  context = g_markup_parse_context_new (&parser, 0, &loader, NULL);

  if (!g_markup_parse_context_parse (context, contents, len, &error))
    {
      g_warning ("Failed to parse ui description: %s", error->message);
      g_markup_parse_context_free (context);
    }

  g_queue_clear (&loader.state);

  free_ux (data);

  data->items = loader.items;
  data->selected_item = loader.items->data;

  data->items_next_id = 0;

  for (l = data->items; l; l = l->next)
    {
      Item *item = l->data;

      if (item->id >= data->items_next_id)
        data->items_next_id = item->id + 1;

      rig_graphable_add_child (data->device_transform, item);
    }

  for (l = loader.entities; l; l = l->next)
    {
      if (rig_graphable_get_parent (l->data) == NULL)
        rig_graphable_add_child (data->scene, l->data);
    }

  g_list_free (data->lights);
  data->lights = loader.lights;

  data->transitions = loader.transitions;
  data->selected_transition = loader.transitions->data;

  data->assets = loader.assets;

  update_asset_list (data);

  rig_shell_queue_redraw (data->ctx->shell);

  g_hash_table_destroy (loader.id_map);

  g_print ("File Loaded\n");
}

static void
init_types (void)
{
  _asset_type_init ();
  _item_type_init ();
  _transition_type_init ();
  _diamond_slice_init_type ();
}

#ifdef __ANDROID__

void
android_main (struct android_app *application)
{
  Data data;

  /* Make sure glue isn't stripped */
  app_dummy ();

  g_android_init ();

  memset (&data, 0, sizeof (Data));
  data.app = application;

  init_types ();

  data.shell = rig_android_shell_new (application,
                                      test_init,
                                      test_fini,
                                      test_paint,
                                      &data);

  data.ctx = rig_context_new (data.shell);

  rig_context_init (data.ctx);

  rig_shell_set_input_callback (data.shell, test_input_handler, &data);

  rig_shell_main (data.shell);
}

#else

int
main (int argc, char **argv)
{
  Data data;
  GOptionContext *context = g_option_context_new (NULL);
  GError *error = NULL;

  g_option_context_add_main_entries (context, rig_handset_entries, NULL);

  if (!g_option_context_parse (context, &argc, &argv, &error))
    {
      g_error ("option parsing failed: %s\n", error->message);
      exit(EXIT_FAILURE);
    }

  memset (&data, 0, sizeof (Data));

  init_types ();

  data.shell = rig_shell_new (test_init, test_fini, test_paint, &data);

  data.ctx = rig_context_new (data.shell);

  rig_context_init (data.ctx);

  rig_shell_set_input_callback (data.shell, test_input_handler, &data);

  rig_shell_main (data.shell);

  return 0;
}

#endif