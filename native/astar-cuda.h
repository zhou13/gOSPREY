#pragma once

#include "astar.h"

#ifdef __cplusplus
extern "C" {
#endif

void init_cuda(void);
int *astar_cuda(bool first_run);

#ifdef __cplusplus
}
#endif
