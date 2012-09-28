
#ifndef _RUT_ENTRY_H_
#define _RUT_ENTRY_H_

#include "rut.h"

extern RutType rut_entry_type;
typedef struct _RutEntry RutEntry;
#define RUT_ENTRY(x) ((RutEntry *) x)

RutEntry *
rut_entry_new (RutContext *context);

void
rut_entry_set_size (RutEntry *entry,
                    float width,
                    float height);

void
rut_entry_set_width (RutEntry *entry,
                     float width);

void
rut_entry_set_height (RutEntry *entry,
                      float height);

void
rut_entry_get_size (RutEntry *entry,
                    float *width,
                    float *height);

RutText *
rut_entry_get_text (RutEntry *entry);

#endif /* _RUT_ENTRY_H_ */