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

#ifndef _RUT_UNDO_JOURNAL_H_
#define _RUT_UNDO_JOURNAL_H_

#include <glib.h>

#include <rut.h>

#include "rig-data.h"
#include "rut-list.h"

typedef enum _UndoRedoOp
{
  UNDO_REDO_PROPERTY_CHANGE_OP,
  UNDO_REDO_N_OPS
} UndoRedoOp;

typedef struct _UndoRedoPropertyChange
{
  RutEntity *entity;
  RutProperty *property;
  RutBoxed value0;
  RutBoxed value1;
} UndoRedoPropertyChange;

typedef struct _UndoRedo
{
  RutList list_node;

  UndoRedoOp op;
  CoglBool mergable;
  union
    {
      UndoRedoPropertyChange prop_change;
    } d;
} UndoRedo;

struct _RigUndoJournal
{
  RigData *data;

  /* List of operations that can be undone. The operations are
   * appended to the end of this list so that they are kept in order
   * from the earliest added operation to the last added operation.
   * The operations are not stored inverted so each operation
   * represents the action that user made. */
  RutList undo_ops;

  /* List of operations that can be redone. As the user presses undo,
   * the operations are added to the tail of this list. Therefore the
   * list is in the order of earliest undone operation to the latest.
   * The operations represent the original action that the user made
   * so it will not need to be inverted before redoing the
   * operation. */
  RutList redo_ops;
};

void
rig_undo_journal_copy_property_and_log (RigUndoJournal *journal,
                                        CoglBool mergable,
                                        RutEntity *entity,
                                        RutProperty *target_prop,
                                        RutProperty *source_prop);

void
rig_undo_journal_log_move (RigUndoJournal *journal,
                           CoglBool mergable,
                           RutEntity *entity,
                           float prev_x,
                           float prev_y,
                           float prev_z,
                           float x,
                           float y,
                           float z);

CoglBool
rig_undo_journal_undo (RigUndoJournal *journal);

CoglBool
rig_undo_journal_redo (RigUndoJournal *journal);

RigUndoJournal *
rig_undo_journal_new (RigData *data);

#endif /* _RUT_UNDO_JOURNAL_H_ */