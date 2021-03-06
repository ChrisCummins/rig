/*
 * Rut.
 *
 * Copyright (C) 2008,2012  Intel Corporation.
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library. If not, see <http://www.gnu.org/licenses/>.
 *
 * Authors:
 *   Øyvind Kolås <pippin@o-hand.com>
 *   Emmanuele Bassi <ebassi@linux.intel.com>
 *   Robert Bragg <robert@linux.intel.com>
 */

#ifndef __RUT_TEXT_H__
#define __RUT_TEXT_H__

#include "rut-types.h"
#include "rut-context.h"
#include "rut-interfaces.h"
#include "rut-text-buffer.h"
#include "rut-closure.h"

#include <pango/pango.h>

G_BEGIN_DECLS

/**
 * SECTION:rut-text
 * @short_description: An actor for displaying and editing text
 *
 * #RutText is an actor that displays custom text using Pango
 * as the text rendering engine.
 *
 * #RutText also allows inline editing of the text if the
 * actor is set editable using rut_text_set_editable().
 *
 * Selection using keyboard or pointers can be enabled using
 * rut_text_set_selectable().
 */


typedef struct _RutText RutText;
#define RUT_TEXT(X) ((RutText *)X)
extern RutType rut_text_type;

RutContext *
rut_text_get_context (RutText *text);

RutTextDirection
rut_text_get_direction (RutText *text);

void
rut_text_set_direction (RutText *text,
                        RutTextDirection direction);

void
_rut_text_init_type (void);

#if 0
void
rut_text_allocate (RutText *text,
                   const RutBox *box);

void
rut_text_get_allocation_box (RutText *text,
                             RutBox *box);
#endif

CoglBool
rut_text_has_overlaps (RutText *text);

typedef void (* RutTextChangedCallback) (RutText *text,
                                         void *user_data);

/**
 * RutText::text-changed:
 * @text: the #RutText that emitted the signal
 *
 * The ::text-changed signal is emitted after @actor's text changes
 */
RutClosure *
rut_text_add_text_changed_callback (RutText *text,
                                    RutTextChangedCallback callback,
                                    void *user_data,
                                    RutClosureDestroyCallback destroy_cb);

typedef void (* RutTextActivateCallback) (RutText *text,
                                          void *user_data);

/**
 * RutText::activate
 * @text: the #RutText that emitted the signal
 *
 * The ::activate signal is emitted each time the actor is 'activated'
 * by the user, normally by pressing the 'Enter' key. The signal is
 * emitted only if #RutText:activatable is set to %TRUE.
 */
RutClosure *
rut_text_add_activate_callback (RutText *text,
                                RutTextActivateCallback callback,
                                void *user_data,
                                RutClosureDestroyCallback destroy_cb);

typedef void (* RutTextCursorEventCallback) (RutText *text,
                                             const RutRectangleInt *rectangle,
                                             void *user_data);
/**
 * RutText::cursor-event:
 * @text: the #RutText that emitted the signal
 * @rectangle: the coordinates of the cursor
 *
 * The ::cursor-event signal is emitted whenever the cursor position
 * changes inside a #RutText actor. Inside @rectangle it is stored
 * the current position and size of the cursor, relative to the actor
 * ittext.
 */
RutClosure *
rut_text_add_cursor_event_callback (RutText *text,
                                    RutTextCursorEventCallback callback,
                                    void *user_data,
                                    RutClosureDestroyCallback destroy_cb);

typedef void (* RutTextInsertedCallback) (RutText *text,
                                          const char *text_str,
                                          int new_text_length,
                                          int *position,
                                          void *user_data);

/**
 * RutText::insert-text:
 * @text: the #RutText that emitted the signal
 * @new_text: the new text to insert
 * @new_text_length: the length of the new text, in bytes, or -1 if
 *     new_text is nul-terminated
 * @position: the position, in characters, at which to insert the
 *     new text. this is an in-out parameter.  After the signal
 *     emission is finished, it should point after the newly
 *     inserted text.
 *
 * This signal is emitted when text is inserted into the actor by
 * the user. It is emitted before @text text changes.
 */
RutClosure *
rut_text_add_text_inserted_callback (RutText *text,
                                     RutTextInsertedCallback callback,
                                     void *user_data,
                                     RutClosureDestroyCallback destroy_cb);

typedef void (* RutTextDeletedCallback) (RutText *text,
                                         int start_pos,
                                         int end_pos,
                                         void *user_data);

/**
 * RutText::delete-text:
 * @text: the #RutText that emitted the signal
 * @start_pos: the starting position
 * @end_pos: the end position
 *
 * This signal is emitted when text is deleted from the actor by
 * the user. It is emitted before @text text changes.
 */
RutClosure *
rut_text_add_text_deleted_callback (RutText *text,
                                    RutTextDeletedCallback callback,
                                    void *user_data,
                                    RutClosureDestroyCallback destroy_cb);

/**
 * rut_text_new:
 *
 * Creates a new #RutText actor. This actor can be used to
 * display and edit text.
 *
 * Return value: the newly created #RutText actor
 */
RutText *
rut_text_new (RutContext *ctx);

/**
 * rut_text_new_full:
 * @font_name: a string with a font description
 * @text: the contents of the actor
 * @buffer: The #RutTextBuffer to use for the text
 *
 * Creates a new #RutText actor, using @font_name as the font
 * description; @text will be used to set the contents of the actor;
 * and @color will be used as the color to render @text.
 *
 * This function is equivalent to calling rut_text_new(),
 * rut_text_set_font_name(), rut_text_set_text() and
 * rut_text_set_color().
 *
 * Return value: the newly created #RutText actor
 */
RutText *
rut_text_new_full (RutContext *ctx,
                   const char *font_name,
                   const char *text,
                   RutTextBuffer *buffer);

/**
 * rut_text_new_with_text:
 * @font_name: (allow-none): a string with a font description
 * @text: the contents of the actor
 *
 * Creates a new #RutText actor, using @font_name as the font
 * description; @text will be used to set the contents of the actor.
 *
 * This function is equivalent to calling rut_text_new(),
 * rut_text_set_font_name(), and rut_text_set_text().
 *
 * Return value: the newly created #RutText actor
 */
RutText *
rut_text_new_with_text (RutContext *ctx,
                        const char *font_name,
                        const char *text);

/**
 * rut_text_new_with_buffer:
 * @buffer: The buffer to use for the new #RutText.
 *
 * Creates a new entry with the specified text buffer.
 *
 * Return value: a new #RutText
 */
RutText *
rut_text_new_with_buffer (RutContext *ctx,
                          RutTextBuffer *buffer);

/**
 * rut_text_get_buffer:
 * @text: a #RutText
 *
 * Get the #RutTextBuffer object which holds the text for
 * this widget.
 *
 * Returns: (transfer none): A #RutTextBuffer object.
 */
RutObject *
rut_text_get_buffer (RutObject *text);

/**
 * rut_text_set_buffer:
 * @text: a #RutText
 * @buffer: a #RutTextBuffer
 *
 * Set the #RutTextBuffer object which holds the text for
 * this widget.
 */
void
rut_text_set_buffer (RutObject *text,
                     RutObject *buffer);

/**
 * rut_text_get_text:
 * @text: a #RutText
 *
 * Retrieves a pointer to the current contents of a #RutText
 * actor.
 *
 * If you need a copy of the contents for manipulating, either
 * use g_strdup() on the returned string, or use:
 *
 * |[
 *    copy = rut_text_get_chars (text, 0, -1);
 * ]|
 *
 * Which will return a newly allocated string.
 *
 * If the #RutText actor is empty, this function will return
 * an empty string, and not %NULL.
 *
 * Return value: (transfer none): the contents of the actor. The returned
 *   string is owned by the #RutText actor and should never be modified
 *   or freed
 */
const char *
rut_text_get_text (RutObject *text);

/**
 * rut_text_set_text:
 * @text: a #RutText
 * @text_str: (allow-none): the text to set. Passing %NULL is the same
 *   as passing "" (the empty string)
 *
 * Sets the contents of a #RutText actor.
 *
 * If the #RutText:use-markup property was set to %TRUE it
 * will be reset to %FALSE as a side effect. If you want to
 * maintain the #RutText:use-markup you should use the
 * rut_text_set_markup() function instead
 */
void
rut_text_set_text (RutObject *text,
                   const char *text_str);

/**
 * rut_text_get_hint_text:
 * @text: a #RutText
 *
 * Retrieves the hint text that will be displayed in the entry box
 * when the text is empty and the entry does not have keyboard focus.
 *
 * Return value: (transfer none): the hint text. The returned
 *   string is owned by the #RutText actor and should never be modified
 *   or freed
 */
const char *
rut_text_get_hint_text (RutObject *text);

/**
 * rut_text_set_hint_text:
 * @text: a #RutText
 * @hint_str: (allow-none): the text to set. Passing %NULL is the same
 *   as passing "" (the empty string)
 *
 * Sets the hint text of a #RutText actor.
 *
 * The hint text is shown when the entry box is empty and it does not
 * have a keyboard focus. It can be used to give a hint to the user
 * about what should be typed in the box.
 */
void
rut_text_set_hint_text (RutObject *text,
                        const char *hint_str);

/**
 * rut_text_set_markup:
 * @text: a #RutText
 * @markup: (allow-none): a string containing Pango markup.
 *   Passing %NULL is the same as passing "" (the empty string)
 *
 * Sets @markup as the contents of a #RutText.
 *
 * This is a convenience function for setting a string containing
 * Pango markup, and it is logically equivalent to:
 *
 * |[
 *   /&ast; the order is important &ast;/
 *   rut_text_set_text (RUT_TEXT (actor), markup);
 *   rut_text_set_use_markup (RUT_TEXT (actor), TRUE);
 * ]|
 */
void
rut_text_set_markup (RutText *text,
                     const char  *markup);

/**
 * rut_text_set_color:
 * @text: a #RutText
 * @color: a #CoglColor
 *
 * Sets the color of the contents of a #RutText actor.
 *
 * The overall opacity of the #RutText actor will be the
 * result of the alpha value of @color and the composited
 * opacity of the actor ittext on the scenegraph, as returned
 * by rut_actor_get_paint_opacity().
 */
void
rut_text_set_color (RutObject *text,
                    const CoglColor *color);

void
rut_text_set_color_u32 (RutText *text,
                        uint32_t u32);

/**
 * rut_text_get_color:
 * @text: a #RutText
 * @color: (out caller-allocates): return location for a #CoglColor
 *
 * Retrieves the text color as set by rut_text_set_color().
 */
const CoglColor *
rut_text_get_color (RutObject *text);

/**
 * rut_text_set_font_name:
 * @text: a #RutText
 * @font_name: (allow-none): a font name, or %NULL to set the default font name
 *
 * Sets the font used by a #RutText. The @font_name string
 * must either be %NULL, which means that the font name from the
 * default #RutBackend will be used; or be something that can
 * be parsed by the pango_font_description_from_string() function,
 * like:
 *
 * |[
 *   rut_text_set_font_name (text, "Sans 10pt");
 *   rut_text_set_font_name (text, "Serif 16px");
 *   rut_text_set_font_name (text, "Helvetica 10");
 * ]|
 */
void
rut_text_set_font_name (RutObject *text,
                        const char *font_name);

/**
 * rut_text_get_font_name:
 * @text: a #RutText
 *
 * Retrieves the font name as set by rut_text_set_font_name().
 *
 * Return value: a string containing the font name. The returned
 *   string is owned by the #RutText actor and should not be
 *   modified or freed
 */
const char *
rut_text_get_font_name (RutObject *text);

/**
 * rut_text_set_font_description:
 * @text: a #RutText
 * @font_desc: a #PangoFontDescription
 *
 * Sets @font_desc as the font description for a #RutText
 *
 * The #PangoFontDescription is copied by the #RutText actor
 * so you can safely call pango_font_description_free() on it after
 * calling this function.
 */
void
rut_text_set_font_description (RutText *text,
                               PangoFontDescription *font_desc);

/**
 * rut_text_get_font_description:
 * @text: a #RutText
 *
 * Retrieves the #PangoFontDescription used by @text
 *
 * Return value: a #PangoFontDescription. The returned value is owned
 *   by the #RutText actor and it should not be modified or freed
 */
PangoFontDescription *
rut_text_get_font_description (RutText *text);

/**
 * rut_text_set_ellipsize:
 * @text: a #RutText
 * @mode: a #PangoEllipsizeMode
 *
 * Sets the mode used to ellipsize (add an ellipsis: "...") to the
 * text if there is not enough space to render the entire contents
 * of a #RutText actor
 */
void
rut_text_set_ellipsize (RutText *text,
                        PangoEllipsizeMode mode);

/**
 * rut_text_get_ellipsize:
 * @text: a #RutText
 *
 * Returns the ellipsizing position of a #RutText actor, as
 * set by rut_text_set_ellipsize().
 *
 * Return value: #PangoEllipsizeMode
 */
PangoEllipsizeMode
rut_text_get_ellipsize (RutText *text);

/**
 * rut_text_set_line_wrap:
 * @text: a #RutText
 * @line_wrap: whether the contents should wrap
 *
 * Sets whether the contents of a #RutText actor should wrap,
 * if they don't fit the size assigned to the actor.
 */
void
rut_text_set_line_wrap (RutObject *text,
                        CoglBool line_wrap);

/**
 * rut_text_get_line_wrap:
 * @text: a #RutText
 *
 * Retrieves the value set using rut_text_set_line_wrap().
 *
 * Return value: %TRUE if the #RutText actor should wrap
 *   its contents
 */
CoglBool
rut_text_get_line_wrap (RutObject *text);

/**
 * rut_text_set_line_wrap_mode:
 * @text: a #RutText
 * @wrap_mode: the line wrapping mode
 *
 * If line wrapping is enabled (see rut_text_set_line_wrap()) this
 * function controls how the line wrapping is performed. The default is
 * %PANGO_WRAP_WORD which means wrap on word boundaries.
 */
void
rut_text_set_line_wrap_mode (RutText *text,
                             PangoWrapMode wrap_mode);

/**
 * rut_text_get_line_wrap_mode:
 * @text: a #RutText
 *
 * Retrieves the line wrap mode used by the #RutText actor.
 *
 * See rut_text_set_line_wrap_mode ().
 *
 * Return value: the wrap mode used by the #RutText
 */
PangoWrapMode
rut_text_get_line_wrap_mode (RutText *text);

/**
 * rut_text_get_layout:
 * @text: a #RutText
 *
 * Retrieves the current #PangoLayout used by a #RutText actor.
 *
 * Return value: (transfer none): a #PangoLayout. The returned object is owned by
 *   the #RutText actor and should not be modified or freed
 */
PangoLayout *
rut_text_get_layout (RutText *text);

/**
 * rut_text_set_attributes:
 * @text: a #RutText
 * @attrs: a #PangoAttrList or %NULL to unset the attributes
 *
 * Sets the attributes list that are going to be applied to the
 * #RutText contents.
 *
 * The #RutText actor will take a reference on the #PangoAttrList
 * passed to this function.
 */
void
rut_text_set_attributes (RutText *text,
                         PangoAttrList *attrs);

/**
 * rut_text_get_attributes:
 * @text: a #RutText
 *
 * Gets the attribute list that was set on the #RutText actor
 * rut_text_set_attributes(), if any.
 *
 * Return value: (transfer none): the attribute list, or %NULL if none was set. The
 *  returned value is owned by the #RutText and should not be unreferenced.
 */
PangoAttrList *
rut_text_get_attributes (RutText *text);

/**
 * rut_text_set_use_markup:
 * @text: a #RutText
 * @setting: %TRUE if the text should be parsed for markup.
 *
 * Sets whether the contents of the #RutText actor contains markup
 * in <link linkend="PangoMarkupFormat">Pango's text markup language</link>.
 *
 * Setting #RutText:use-markup on an editable #RutText will
 * not have any effect except hiding the markup.
 *
 * See also #RutText:use-markup.
 */
void
rut_text_set_use_markup (RutObject *text,
                         CoglBool setting);

/**
 * rut_text_get_use_markup:
 * @text: a #RutText
 *
 * Retrieves whether the contents of the #RutText actor should be
 * parsed for the Pango text markup.
 *
 * Return value: %TRUE if the contents will be parsed for markup
 */
CoglBool
rut_text_get_use_markup (RutObject *text);

/**
 * rut_text_set_line_alignment:
 * @text: a #RutText
 * @alignment: A #PangoAlignment
 *
 * Sets the way that the lines of a wrapped label are aligned with
 * respect to each other. This does not affect the overall alignment
 * of the label within its allocated or specified width.
 *
 * To align a #RutText actor you should add it to a container
 * that supports alignment, or use the anchor point.
 */
void
rut_text_set_line_alignment (RutText *text,
                             PangoAlignment alignment);

/**
 * rut_text_get_line_alignment:
 * @text: a #RutText
 *
 * Retrieves the alignment of a #RutText, as set by
 * rut_text_set_line_alignment().
 *
 * Return value: a #PangoAlignment
 */
PangoAlignment
rut_text_get_line_alignment (RutText *text);

/**
 * rut_text_set_justify:
 * @text: a #RutText
 * @justify: whether the text should be justified
 *
 * Sets whether the text of the #RutText actor should be justified
 * on both margins. This setting is ignored if Rut is compiled
 * against Pango &lt; 1.18.
 */
void
rut_text_set_justify (RutObject *text,
                      CoglBool justify);

/**
 * rut_text_get_justify:
 * @text: a #RutText
 *
 * Retrieves whether the #RutText actor should justify its contents
 * on both margins.
 *
 * Return value: %TRUE if the text should be justified
 */
CoglBool
rut_text_get_justify (RutObject *text);

/**
 * rut_text_insert_unichar:
 * @text: a #RutText
 * @wc: a Unicode character
 *
 * Inserts @wc at the current cursor position of a
 * #RutText actor.
 */
void
rut_text_insert_unichar (RutText *text,
                         uint32_t wc);

/**
 * rut_text_delete_chars:
 * @text: a #RutText
 * @n_chars: the number of characters to delete
 *
 * Deletes @n_chars inside a #RutText actor, starting from the
 * current cursor position.
 *
 * Somewhat awkwardly, the cursor position is decremented by the same
 * number of characters you've deleted.
 */
void
rut_text_delete_chars (RutText *text,
                       unsigned int n_chars);

/**
 * rut_text_insert_text:
 * @text: a #RutText
 * @text: the text to be inserted
 * @position: the position of the insertion, or -1
 *
 * Inserts @text into a #RutActor at the given position.
 *
 * If @position is a negative number, the text will be appended
 * at the end of the current contents of the #RutText.
 *
 * The position is expressed in characters, not in bytes.
 */
void
rut_text_insert_text (RutText *text,
                      const char *text_str,
                      int position);

/**
 * rut_text_delete_text:
 * @text: a #RutText
 * @start_pos: starting position
 * @end_pos: ending position
 *
 * Deletes the text inside a #RutText actor between @start_pos
 * and @end_pos.
 *
 * The starting and ending positions are expressed in characters,
 * not in bytes.
 */
void
rut_text_delete_text (RutText *text,
                      int start_pos,
                      int end_pos);

/**
 * rut_text_get_chars:
 * @text: a #RutText
 * @start_pos: start of text, in characters
 * @end_pos: end of text, in characters
 *
 * Retrieves the contents of the #RutText actor between
 * @start_pos and @end_pos, but not including @end_pos.
 *
 * The positions are specified in characters, not in bytes.
 *
 * Return value: a newly allocated string with the contents of
 *   the text actor between the specified positions. Use g_free()
 *   to free the resources when done
 */
char *
rut_text_get_chars (RutText *text,
                    int start_pos,
                    int end_pos);

/**
 * rut_text_set_editable:
 * @text: a #RutText
 * @editable: whether the #RutText should be editable
 *
 * Sets whether the #RutText actor should be editable.
 *
 * An editable #RutText with key focus set using
 * rut_actor_grab_key_focus() or rut_stage_set_key_focus()
 * will receive key events and will update its contents accordingly.
 */
void
rut_text_set_editable (RutObject *text,
                       CoglBool editable);

/**
 * rut_text_get_editable:
 * @text: a #RutText
 *
 * Retrieves whether a #RutText is editable or not.
 *
 * Return value: %TRUE if the actor is editable
 */
CoglBool
rut_text_get_editable (RutObject *text);

/**
 * rut_text_set_activatable:
 * @text: a #RutText
 * @activatable: whether the #RutText actor should be activatable
 *
 * Sets whether a #RutText actor should be activatable.
 *
 * An activatable #RutText actor will emit the #RutText::activate
 * signal whenever the 'Enter' (or 'Return') key is pressed; if it is not
 * activatable, a new line will be appended to the current content.
 *
 * An activatable #RutText must also be set as editable using
 * rut_text_set_editable().
 */
void
rut_text_set_activatable (RutObject *text,
                          CoglBool activatable);

/**
 * rut_text_get_activatable:
 * @text: a #RutText
 *
 * Retrieves whether a #RutText is activatable or not.
 *
 * Return value: %TRUE if the actor is activatable
 */
CoglBool
rut_text_get_activatable (RutObject *text);

/**
 * rut_text_get_cursor_position:
 * @text: a #RutText
 *
 * Retrieves the cursor position.
 *
 * Return value: the cursor position, in characters
 */
int
rut_text_get_cursor_position (RutObject *text);

/**
 * rut_text_set_cursor_position:
 * @text: a #RutText
 * @position: the new cursor position, in characters
 *
 * Sets the cursor of a #RutText actor at @position.
 *
 * The position is expressed in characters, not in bytes.
 */
void
rut_text_set_cursor_position (RutObject *text,
                              int position);

/**
 * rut_text_set_cursor_visible:
 * @text: a #RutText
 * @cursor_visible: whether the cursor should be visible
 *
 * Sets whether the cursor of a #RutText actor should be
 * visible or not.
 *
 * The color of the cursor will be the same as the text color
 * unless rut_text_set_cursor_color() has been called.
 *
 * The size of the cursor can be set using rut_text_set_cursor_size().
 *
 * The position of the cursor can be changed programmatically using
 * rut_text_set_cursor_position().
 */
void
rut_text_set_cursor_visible (RutObject *text,
                             CoglBool cursor_visible);

/**
 * rut_text_get_cursor_visible:
 * @text: a #RutText
 *
 * Retrieves whether the cursor of a #RutText actor is visible.
 *
 * Return value: %TRUE if the cursor is visible
 */
CoglBool
rut_text_get_cursor_visible (RutObject *text);

/**
 * rut_text_set_cursor_color:
 * @text: a #RutText
 * @color: the color of the cursor, or %NULL to unset it
 *
 * Sets the color of the cursor of a #RutText actor.
 *
 * If @color is %NULL, the cursor color will be the same as the
 * text color.
 */
void
rut_text_set_cursor_color (RutObject *text,
                           const CoglColor *color);

void
rut_text_set_cursor_color_u32 (RutText *text,
                               uint32_t u32);

/**
 * rut_text_get_cursor_color:
 * @text: a #RutText
 * @color: (out): return location for a #CoglColor
 *
 * Retrieves the color of the cursor of a #RutText actor.
 */
const CoglColor *
rut_text_get_cursor_color (RutObject *text);

CoglBool
rut_text_get_cursor_color_set (RutObject *text);

/**
 * rut_text_set_cursor_size:
 * @text: a #RutText
 * @size: the size of the cursor, in pixels, or -1 to use the
 *   default value
 *
 * Sets the size of the cursor of a #RutText. The cursor
 * will only be visible if the #RutText:cursor-visible property
 * is set to %TRUE.
 */
void
rut_text_set_cursor_size (RutObject *text,
                          int size);

/**
 * rut_text_get_cursor_size:
 * @text: a #RutText
 *
 * Retrieves the size of the cursor of a #RutText actor.
 *
 * Return value: the size of the cursor, in pixels
 */
int
rut_text_get_cursor_size (RutObject *text);

/**
 * rut_text_set_selectable:
 * @text: a #RutText
 * @selectable: whether the #RutText actor should be selectable
 *
 * Sets whether a #RutText actor should be selectable.
 *
 * A selectable #RutText will allow selecting its contents using
 * the pointer or the keyboard.
 */
void
rut_text_set_selectable (RutObject *text,
                         CoglBool selectable);

/**
 * rut_text_get_selectable:
 * @text: a #RutText
 *
 * Retrieves whether a #RutText is selectable or not.
 *
 * Return value: %TRUE if the actor is selectable
 */
CoglBool
rut_text_get_selectable (RutObject *text);

/**
 * rut_text_set_selection_bound:
 * @text: a #RutText
 * @selection_bound: the position of the end of the selection, in characters
 *
 * Sets the other end of the selection, starting from the current
 * cursor position.
 *
 * If @selection_bound is -1, the selection unset.
 */
void
rut_text_set_selection_bound (RutObject *text,
                              int selection_bound);

/**
 * rut_text_get_selection_bound:
 * @text: a #RutText
 *
 * Retrieves the other end of the selection of a #RutText actor,
 * in characters from the current cursor position.
 *
 * Return value: the position of the other end of the selection
 */
int
rut_text_get_selection_bound (RutObject *text);

/**
 * rut_text_set_selection:
 * @text: a #RutText
 * @start_pos: start of the selection, in characters
 * @end_pos: end of the selection, in characters
 *
 * Selects the region of text between @start_pos and @end_pos.
 *
 * This function changes the position of the cursor to match
 * @start_pos and the selection bound to match @end_pos.
 */
void
rut_text_set_selection (RutText *text,
                        int start_pos,
                        int end_pos);

/**
 * rut_text_get_selection:
 * @text: a #RutText
 *
 * Retrieves the currently selected text.
 *
 * Return value: a newly allocated string containing the currently
 *   selected text, or %NULL. Use g_free() to free the returned
 *   string.
 */
char *
rut_text_get_selection (RutText *text);

/**
 * rut_text_set_selection_color:
 * @text: a #RutText
 * @color: the color of the selection, or %NULL to unset it
 *
 * Sets the color of the selection of a #RutText actor.
 *
 * If @color is %NULL, the selection color will be the same as the
 * cursor color, or if no cursor color is set either then it will be
 * the same as the text color.
 */
void
rut_text_set_selection_color (RutObject *text,
                              const CoglColor *color);

void
rut_text_set_selection_color_u32 (RutText *text,
                                  uint32_t u32);

/**
 * rut_text_get_selection_color:
 * @text: a #RutText
 * @color: (out caller-allocates): return location for a #CoglColor
 *
 * Retrieves the color of the selection of a #RutText actor.
 */
const CoglColor *
rut_text_get_selection_color (RutObject *text);

CoglBool
rut_text_get_selection_color_set (RutObject *text);

/**
 * rut_text_delete_selection:
 * @text: a #RutText
 *
 * Deletes the currently selected text
 *
 * This function is only useful in subclasses of #RutText
 *
 * Return value: %TRUE if text was deleted or if the text actor
 *   is empty, and %FALSE otherwise
 */
CoglBool
rut_text_delete_selection (RutText *text);

/**
 * rut_text_set_password_char:
 * @text: a #RutText
 * @wc: a Unicode character, or 0 to unset the password character
 *
 * Sets the character to use in place of the actual text in a
 * password text actor.
 *
 * If @wc is 0 the text will be displayed as it is entered in the
 * #RutText actor.
 */
void
rut_text_set_password_char (RutObject *text,
                            uint32_t wc);

/**
 * rut_text_get_password_char:
 * @text: a #RutText
 *
 * Retrieves the character to use in place of the actual text
 * as set by rut_text_set_password_char().
 *
 * Return value: a Unicode character or 0 if the password
 *   character is not set
 */
uint32_t
rut_text_get_password_char (RutObject *text);

/**
 * rut_text_set_max_length:
 * @text: a #RutText
 * @max: the maximum number of characters allowed in the text actor; 0
 *   to disable or -1 to set the length of the current string
 *
 * Sets the maximum allowed length of the contents of the actor. If the
 * current contents are longer than the given length, then they will be
 * truncated to fit.
 */
void
rut_text_set_max_length (RutObject *text,
                         int max);

/**
 * rut_text_get_max_length:
 * @text: a #RutText
 *
 * Gets the maximum length of text that can be set into a text actor.
 *
 * See rut_text_set_max_length().
 *
 * Return value: the maximum number of characters.
 */
int
rut_text_get_max_length (RutObject *text);

/**
 * rut_text_set_single_line_mode:
 * @text: a #RutText
 * @single_line: whether to enable single line mode
 *
 * Sets whether a #RutText actor should be in single line mode
 * or not. Only editable #RutText<!-- -->s can be in single line
 * mode.
 *
 * A text actor in single line mode will not wrap text and will clip
 * the visible area to the predefined size. The contents of the
 * text actor will scroll to display the end of the text if its length
 * is bigger than the allocated width.
 *
 * When setting the single line mode the #RutText:activatable
 * property is also set as a side effect. Instead of entering a new
 * line character, the text actor will emit the #RutText::activate
 * signal.
 */
void
rut_text_set_single_line_mode (RutObject *text,
                               CoglBool single_line);

/**
 * rut_text_get_single_line_mode:
 * @text: a #RutText
 *
 * Retrieves whether the #RutText actor is in single line mode.
 *
 * Return value: %TRUE if the #RutText actor is in single line mode
 */
CoglBool
rut_text_get_single_line_mode (RutObject *text);

/**
 * rut_text_set_selected_text_color:
 * @text: a #RutText
 * @color: the selected text color, or %NULL to unset it
 *
 * Sets the selected text color of a #RutText actor.
 *
 * If @color is %NULL, the selected text color will be the same as the
 * selection color, which then falls back to cursor, and then text color.
 */
void
rut_text_set_selected_text_color (RutObject *text,
                                  const CoglColor *color);

void
rut_text_set_selected_text_color_u32 (RutText *text,
                                      uint32_t u32);

/**
 * rut_text_get_selected_text_color:
 * @text: a #RutText
 * @color: (out caller-allocates): return location for a #CoglColor
 *
 * Retrieves the color of selected text of a #RutText actor.
 */
const CoglColor *
rut_text_get_selected_text_color (RutObject *text);

CoglBool
rut_text_get_selected_text_color_set (RutObject *text);

/**
 * rut_text_activate:
 * @text: a #RutText
 *
 * Emits the #RutText::activate signal, if @text has been set
 * as activatable using rut_text_set_activatable().
 *
 * This function can be used to emit the ::activate signal inside
 * a #RutActor::captured-event or #RutActor::key-press-event
 * signal handlers before the default signal handler for the
 * #RutText is invoked.
 *
 * Return value: %TRUE if the ::activate signal has been emitted,
 *   and %FALSE otherwise
 */
CoglBool
rut_text_activate (RutText *text);

/**
 * rut_text_coords_to_position:
 * @text: a #RutText
 * @x: the X coordinate, relative to the actor
 * @y: the Y coordinate, relative to the actor
 *
 * Retrieves the position of the character at the given coordinates.
 *
 * Return: the position of the character
 */
int
rut_text_coords_to_position (RutText *text,
                             float x,
                             float y);

/**
 * rut_text_position_to_coords:
 * @text: a #RutText
 * @position: position in characters
 * @x: (out): return location for the X coordinate, or %NULL
 * @y: (out): return location for the Y coordinate, or %NULL
 * @line_height: (out): return location for the line height, or %NULL
 *
 * Retrieves the coordinates of the given @position.
 *
 * Return value: %TRUE if the conversion was successful
 */
CoglBool
rut_text_position_to_coords (RutText *text,
                             int position,
                             float *x,
                             float *y,
                             float *line_height);

/**
 * rut_text_set_preedit_string:
 * @text: a #RutText
 * @preedit_str: (allow-none): the pre-edit string, or %NULL to unset it
 * @preedit_attrs: (allow-none): the pre-edit string attributes
 * @cursor_pos: the cursor position for the pre-edit string
 *
 * Sets, or unsets, the pre-edit string. This function is useful
 * for input methods to display a string (with eventual specific
 * Pango attributes) before it is entered inside the #RutText
 * buffer.
 *
 * The preedit string and attributes are ignored if the #RutText
 * actor is not editable.
 *
 * This function should not be used by applications
 */
void
rut_text_set_preedit_string (RutText *text,
                             const char *preedit_str,
                             PangoAttrList *preedit_attrs,
                             unsigned int cursor_pos);

/**
 * rut_text_get_layout_offsets:
 * @text: a #RutText
 * @x: (out): location to store X offset of layout, or %NULL
 * @y: (out): location to store Y offset of layout, or %NULL
 *
 * Obtains the coordinates where the #RutText will draw the #PangoLayout
 * representing the text.
 */
void
rut_text_get_layout_offsets (RutText *text,
                             int *x,
                             int *y);

/**
 * rut_text_grab_key_focus:
 * @text: a #RutText
 *
 * Causes @text to try to grab the key focus.
 */
void
rut_text_grab_key_focus (RutText *text);

void
rut_text_ungrab_key_focus (RutText *text);

void
rut_text_set_width (RutObject *text, float width);

void
rut_text_set_height (RutObject *text, float height);

RutMesh *
rut_text_get_pick_mesh (RutText *text);

G_END_DECLS

#endif /* __RUT_TEXT_H__ */
