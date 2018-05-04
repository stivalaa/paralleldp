/*****************************************************************************
 * 
 * File:    cellpool.c
 * Author:  Alex Stivala
 * Created: April 2009
 *
 * Simple cell pool allocator. By allocating a big chunk of memory at
 * initialization time, we can just get one cell at a time from it
 * in a lock-free manner without having to call malloc() again,
 * so that multiple threads can allocate cells from the pool without
 * locking.
 *
 * On Linux,
 * gcc version 4.1.0 or greater is required, in order to use the
 * __sync_val_compare_and_swap() builtin
 * (this module was developed on Linux 2.6.22 (x86) with gcc 4.1.3).
 * On Solaris, we use the umem_cache_*() functions instead, rather
 * than implement a cell pool allocator ourselves.
 * (Actually, we don't since it turns out to be slower+less scalable than
 * the trivial cellpool implemented here). Define USE_SOLARIS_UMEM to use
 * the umem functions on Solaris.
 *
 * TODO don't even have  facility to free/reuse cells, don't need it
 * since we just terminate the process.
 * TODO only one cell pool allowed (file static data). Should change
 * to have handles so user can allocate different pools.
 * 
 *
 * $Id: cellpool.c 2324 2009-05-06 02:44:04Z astivala $
 *
 *****************************************************************************/

#undef USE_SOLARIS_UMEM /* slower than trivial celpool,. don't use */

#include <stdlib.h>
#include "cellpool.h"
#ifdef USE_SOLARIS_UMEM
#include <umem.h>
#else
#include "atomicdefs.h"
#endif

#ifdef USE_SOLARIS_UMEM
static umem_cache_t *cache; /* The Solaris umem cache */
#else
static void  *cellpool; /* The cell pool. */
static size_t cellsize; /* size of a cell */
static size_t poolsize; /* total size of cell pool */
static void  *nextcell; /* pointer to next cell in the pool
                          serialized with CAS logic in cellpool_alloc() */
#endif


/*
 * cellpool_alloc()
 *
 *   allocate a new cell of previously intialized size from pool
 *
 *   Parameters:
 *      None.
 *
 *   Return value:
 *      Pointer to new cell, or NULL if none available.
 *
 *   Uses static data:
 *      nextcell -pointer to next cell in cellpool (read/write, CAS serialized)
 */
void *cellpool_alloc(void)
{
#ifdef USE_SOLARIS_UMEM
  return umem_cache_alloc(cache, UMEM_DEFAULT);
#else
  void *cell, *newnextcell;
#ifdef USE_THREADING
  do
  {
    if ((char *)nextcell >= (char *)cellpool + poolsize)
      return NULL;
    cell = nextcell;
    newnextcell = (char *)cell + cellsize;
  }
  while (CASPTR(&nextcell,cell,newnextcell) != cell);
#else
  if ((char *)nextcell >= (char *)cellpool + poolsize)
    return NULL;
  cell = nextcell;
  newnextcell = (char *)cell + cellsize;
  nextcell = newnextcell;
#endif
  return cell;
#endif 
}

/*
 * cellpool_initialize()
 *
 *   initialize the cell pool
 *
 *   Parameters:
 *     cell_size - Size of each cell
 *     num_cells - Number of cells to allocate
 *
 *   Return value:
 *     Pointer to start of cell pool memory (first cell) or NULL
 *     on failure
 *
 *    Uses static data:
 *       cellpool - the cell pool (write)          
 *       cellsize - size of each cell (write)
 *       poolsize - total size of cell pool (write)
 *       nextcell - pointer to next cell in cellpool (write)
 */
void *cellpool_initialize(size_t cell_size, int num_cells)
{
#ifdef USE_SOLARIS_UMEM
  cache = umem_cache_create("cellpool", cell_size, 0, 
                            NULL, NULL, NULL, NULL, NULL, 0);
  return cache;
#else
  cellsize = cell_size;
  poolsize = num_cells * cell_size;
  if (!(cellpool = malloc(poolsize)))
    return NULL;
  nextcell = cellpool;
  return cellpool;
#endif
}
