#ifndef _RUT_TIMELINE_H_
#define _RUT_TIMELINE_H_

void
_rut_timeline_init_type (void);

typedef struct _RutTimeline RutTimeline;
#define RUT_TIMELINE(x) ((RutTimeline *)x)

RutTimeline *
rut_timeline_new (RutContext *ctx,
                  float length);

void
rut_timeline_start (RutTimeline *timeline);

void
rut_timeline_stop (RutTimeline *timeline);

CoglBool
rut_timeline_is_running (RutTimeline *timeline);

double
rut_timeline_get_elapsed (RutObject *timeline);

void
rut_timeline_set_elapsed (RutObject *timeline,
                          double elapsed);

double
rut_timeline_get_progress (RutObject *timeline);

void
rut_timeline_set_progress (RutObject *timeline,
                           double elapsed);

void
rut_timeline_set_loop_enabled (RutTimeline *timeline, CoglBool enabled);

/* PRIVATE */
void
_rut_timeline_update (RutTimeline *timeline);

#endif /* _RUT_TIMELINE_H_ */
