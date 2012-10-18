/*
 * Rut
 *
 * A tiny toolkit
 *
 * Copyright (C) 2012 Intel Corporation.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
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

#ifndef _RUT_UI_VIEWPORT_H_
#define _RUT_UI_VIEWPORT_H_

#include <cogl/cogl.h>

#include "rut-type.h"

typedef struct _RutUIViewport RutUIViewport;
#define RUT_UI_VIEWPORT(X) ((RutUIViewport *)X)

extern RutType rut_ui_viewport_type;

RutUIViewport *
rut_ui_viewport_new (RutContext *ctx,
                     float width,
                     float height,
                     ...);

void
rut_ui_viewport_set_width (RutUIViewport *ui_viewport, float width);

void
rut_ui_viewport_set_height (RutUIViewport *ui_viewport, float height);

void
rut_ui_viewport_set_size (RutUIViewport *ui_viewport,
                          float width,
                          float height);

void
rut_ui_viewport_get_size (RutUIViewport *ui_viewport,
                          float *width,
                          float *height);

void
rut_ui_viewport_set_doc_x (RutUIViewport *ui_viewport, float doc_x);

void
rut_ui_viewport_set_doc_y (RutUIViewport *ui_viewport, float doc_y);

void
rut_ui_viewport_set_doc_width (RutUIViewport *ui_viewport, float doc_width);

void
rut_ui_viewport_set_doc_height (RutUIViewport *ui_viewport, float doc_height);

void
rut_ui_viewport_set_doc_scale_x (RutUIViewport *ui_viewport, float doc_scale_x);

void
rut_ui_viewport_set_doc_scale_y (RutUIViewport *ui_viewport, float doc_scale_y);

float
rut_ui_viewport_get_width (RutUIViewport *ui_viewport);

float
rut_ui_viewport_get_height (RutUIViewport *ui_viewport);

float
rut_ui_viewport_get_doc_x (RutUIViewport *ui_viewport);

float
rut_ui_viewport_get_doc_y (RutUIViewport *ui_viewport);

float
rut_ui_viewport_get_doc_scale_x (RutUIViewport *ui_viewport);

float
rut_ui_viewport_get_doc_scale_y (RutUIViewport *ui_viewport);

const CoglMatrix *
rut_ui_viewport_get_doc_matrix (RutUIViewport *ui_viewport);

RutObject *
rut_ui_viewport_get_doc_node (RutUIViewport *ui_viewport);

void
rut_ui_viewport_set_x_pannable (RutUIViewport *ui_viewport,
                                CoglBool pannable);

CoglBool
rut_ui_viewport_get_x_pannable (RutUIViewport *ui_viewport);

void
rut_ui_viewport_set_y_pannable (RutUIViewport *ui_viewport,
                                CoglBool pannable);

CoglBool
rut_ui_viewport_get_y_pannable (RutUIViewport *ui_viewport);

/**
 * rut_ui_viewport_set_sync_widget:
 * @ui_viewport: The viewport object
 * @widget: An object implementing the RutSizable interface
 *
 * Sets a widget to use to specify the doc size. The viewport will
 * track the preferred size of the widget and set the doc to the same
 * size whenever it changes. It will also then set the widget to its
 * preferred size, or larger if the expand properties are enabled. The
 * widget should typically be a child of the doc.
 */
void
rut_ui_viewport_set_sync_widget (RutUIViewport *ui_viewport,
                                 RutObject *widget);

/**
 * rut_ui_viewport_set_x_expand:
 * @ui_viewport: The viewport object
 * @expand: The new value
 *
 * Sets whether the viewport should set the doc size to be wider than
 * the sync widget's preferred size if the viewport itself is wider.
 * This only has any effect if the sync_widget property is set.
 */
void
rut_ui_viewport_set_x_expand (RutUIViewport *ui_viewport,
                              CoglBool expand);

/**
 * rut_ui_viewport_set_y_expand:
 * @ui_viewport: The viewport object
 * @expand: The new value
 *
 * Sets whether the viewport should set the doc size to be taller than
 * the sync widget's preferred size if the viewport itself is taller.
 * This only has any effect if the sync_doc_size property is enabled.
 */
void
rut_ui_viewport_set_y_expand (RutUIViewport *ui_viewport,
                              CoglBool expand);

#endif /* _RUT_UI_VIEWPORT_H_ */