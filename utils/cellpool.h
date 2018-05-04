#ifndef CELLPOOL_H
#define CELLPOOL_H
/*****************************************************************************
 * 
 * File:    cellpool.h
 * Author:  Alex Stivala
 * Created: April 2009
 *
 * Declarations for simple cell pool allocator.
 * 
 *
 * $Id: cellpool.h 2264 2009-04-22 02:49:57Z astivala $
 *
 *****************************************************************************/

/* allocate a new cell of previously intialized size from pool */
void *cellpool_alloc(void);

/* initizlie the cell pool */
void *cellpool_initialize(size_t cell_size, int num_cells);

#endif /* CELLPOOL_H */
