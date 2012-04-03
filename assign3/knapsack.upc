#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <upc.h>

//#pragma upc strict


//
// Constant Definitions
//
#define CAPACITY    31999
#define NITEMS      16000
//#define CAPACITY    20
//#define NITEMS      4
#define BLK_WIDTH   ((CAPACITY+1+THREADS-1) / THREADS)

//
// Type Definitions
//    
typedef shared [BLK_WIDTH]  int* sintptr_cblk;
typedef shared [NITEMS]     int* sintptr_nblk;
typedef volatile strict shared [1] int* sintptr_sngl;

//
// auxiliary functions
//
inline int max( int a, int b ) { return a > b ? a : b; }
inline int min( int a, int b ) { return a < b ? a : b; }

double read_timer( )
{
    static int initialized = 0;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = 1;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}


//
//  serial solver to check correctness
//
int solve_serial( int nitems, int cap, int *w, int *v, int* T )
{
    int i, j, best, wj, vj;
    
    //build_table locally
    wj = w[0];
    vj = v[0];
    for( i = 0;  i <  wj;  i++ ) T[i] = 0;
    for( i = wj; i <= cap; i++ ) T[i] = vj;
    for( j = 1; j < nitems; j++ ) 
    {
        wj = w[j];
        vj = v[j];
        for( i = 0;  i <  wj;  i++ ) T[i+cap+1] = T[i];
        for( i = wj; i <= cap; i++ ) T[i+cap+1] = max( T[i], T[i-wj]+vj );
        T += cap+1;
    }
    best = T[cap];
    
    return best;
}

void backtrack_serial( int nitems, int cap, int *T, int *w, int *u )
{
    int i, j;
    
    i = nitems*(cap+1) - 1;
    for( j = nitems-1; j > 0; j-- )
    {
        u[j] = T[i] != T[i-cap-1];
        i -= cap+1 + (u[j] ? w[j] : 0 );
    }
    u[0] = T[i] != 0;
}

void compute_serial(int nitems, int cap, int *w, int *v, 
                    int* nused, int* weight, int* value){
    //alloc local resources
    int* T = malloc( nitems*(cap+1)*sizeof(int) );
    int* used = malloc( nitems*sizeof(int) );
    if( !T )
    {
        fprintf( stderr, "Failed to allocate memory" );
        upc_global_exit( -1 );
    }

    *value = solve_serial(nitems, cap, w, v, T);
    backtrack_serial(nitems, cap, T, w, used);

    // Determine total weight and value, and num items used, by iterating through UPC output array
    *weight = 0; 
    *nused = 0;
    for( int i = 0; i < nitems; i++ )
    {
        if( used[i] )
        {
            (*nused)++;
            *weight += w[i];
        }
    }
}


//
// Solves the knapsack problem using shared memory
//
int solve_upc(  int nitems, int cap, sintptr_cblk s_old_row, sintptr_cblk s_new_row, sintptr_cblk s_old_cnt,
                sintptr_cblk s_new_cnt, int* l_old_row, int* l_new_row, int* l_old_cnt, int* l_new_cnt,
                int* scratch, int* scratch_cnt, int* l_weight, int* l_value, int* nused, int* total_weight,
                sintptr_sngl s_progress)
{
    // Local variables
    int i, j;
    int start_col, scratch_start_col, scratch_end_col, boundary_col;
    int *temp_intptr;
    int value_with_item, cnt_with_item;
    sintptr_cblk temp_sintptr;

    // Intialize first row
    int val_item0 = l_value[0];

    upc_forall(i = 0; i < l_weight[0]; i++; &s_old_row[i])
    {
        s_old_row[i] = 0;
        s_old_cnt[i] = 0;
    }

    upc_forall(i = l_weight[0]; i < cap+1; i++; &s_old_row[i])
    {
        s_old_row[i] = val_item0;
        s_old_cnt[i] = 1;
    }

    s_progress[MYTHREAD] = 0;

    upc_barrier;

    // Determine this thread's starting column
    start_col = MYTHREAD * BLK_WIDTH;

    // Iterate through remaining rows
    for(i = 1; i < nitems; i++)
    {
        // Copy the section of the above row, that we need and that belongs to another thread, to scratch
        // Note: before using scratch, took 15 seconds
        scratch_start_col = max(0, start_col - l_weight[i]);
        scratch_end_col = min( max(0, start_col - l_weight[i] + BLK_WIDTH) , start_col);  // end is exclusive
        boundary_col = min( (scratch_start_col / BLK_WIDTH + 1) * BLK_WIDTH, scratch_end_col );

//        for(j=scratch_start_col; j < scratch_end_col; j++)
//            scratch[j-scratch_start_col] = s_old_row[j];

        upc_memget(scratch, &s_old_row[scratch_start_col], (boundary_col - scratch_start_col)*sizeof(int));
        upc_memget(&scratch[boundary_col - scratch_start_col], &s_old_row[boundary_col], (scratch_end_col - boundary_col)*sizeof(int));

        upc_memget(scratch_cnt, &s_old_cnt[scratch_start_col], (boundary_col - scratch_start_col)*sizeof(int));
        upc_memget(&scratch_cnt[boundary_col - scratch_start_col], &s_old_cnt[boundary_col], (scratch_end_col - boundary_col)*sizeof(int));
/*
        upc_memget(scratch, &s_old_row[scratch_start_col], 1*sizeof(int));
        upc_memget(&scratch[boundary_col - scratch_start_col], &s_old_row[boundary_col], 1*sizeof(int));

        upc_memget(scratch_cnt, &s_old_cnt[scratch_start_col], 1*sizeof(int));
        upc_memget(&scratch_cnt[boundary_col - scratch_start_col], &s_old_cnt[boundary_col], 1*sizeof(int));
*/
        // Iterate through all columns belonging to this thread
        // Note: the last thread does some extra work on non-existent columns
        for(j = 0; j < BLK_WIDTH; j++)
        {
            // If this item is larger than the current capacity
            if(start_col + j < l_weight[i])
            {
                l_new_row[j] = l_old_row[j];
                l_new_cnt[j] = l_old_cnt[j];
            }

            else
            {
                if(j - l_weight[i] < 0)
                {
                    value_with_item = scratch[start_col + j - l_weight[i] - scratch_start_col] + l_value[i];
                    cnt_with_item = scratch_cnt[start_col + j - l_weight[i] - scratch_start_col] + 1;
                }
                else
                {
                    value_with_item = l_old_row[j - l_weight[i]] + l_value[i];
                    cnt_with_item = l_old_cnt[j - l_weight[i]] + 1;
                }

                // If not using the item
                if(l_old_row[j] >= value_with_item)
                {
                    l_new_row[j] = l_old_row[j];
                    l_new_cnt[j] = l_old_cnt[j];
                }
                // Else if using the item
                else
                {
                    l_new_row[j] = value_with_item;
                    l_new_cnt[j] = cnt_with_item;
                }
            }
        }

        s_progress[MYTHREAD] = i;

        // Sync all threads
        //upc_barrier;
        if(i == nitems - 1)
            upc_barrier;
        else
        {
            scratch_start_col = max(0, start_col - l_weight[i+1]);
            scratch_end_col = min( max(0, start_col - l_weight[i+1] + BLK_WIDTH) , start_col);  // end is exclusive
            int dep_start_col = start_col + l_weight[i];
            int dep_end_col = (dep_start_col + BLK_WIDTH);  // end is exclusive

            int t1 = scratch_start_col / BLK_WIDTH;
            int t2 = max(0, (scratch_end_col - 1) / BLK_WIDTH);
            int t3 = min(THREADS-1, dep_start_col / BLK_WIDTH);
            int t4 = min(THREADS-1, (dep_end_col - 1) / BLK_WIDTH);
            while(s_progress[t1] < i) ;
            while(s_progress[t2] < i) ;
            while(s_progress[t3] < i) ;
            while(s_progress[t4] < i) ;
        }

        // Swap new and old rows
        temp_intptr = l_old_row;
        l_old_row = l_new_row;
        l_new_row = temp_intptr;

        temp_sintptr = s_old_row;
        s_old_row = s_new_row;
        s_new_row = temp_sintptr;

        temp_intptr = l_old_cnt;
        l_old_cnt = l_new_cnt;
        l_new_cnt = temp_intptr;

        temp_sintptr = s_old_cnt;
        s_old_cnt = s_new_cnt;
        s_new_cnt = temp_sintptr;
    }

    // Calculate total weight
    if(MYTHREAD == THREADS-1)
    {
        *total_weight = -1;

        for(i = CAPACITY; i > 0; i--)
        {
            if(s_old_row[i] != s_old_row[i-1])
            {
                *total_weight = i;
                break;
            }
        }

        // Transfer weight from last thread to first
        s_new_cnt[0] = *total_weight;   
    }

    upc_barrier;

    // Return results
    if(MYTHREAD == 0)
    {
        *nused = s_old_cnt[cap];
        *total_weight = s_new_cnt[0];
        return s_old_row[cap];
    }

    else
        return 0;
}


//
//  benchmarking program
//
int main( int argc, char** argv )
{
    // Local variables
    int i, best_value, best_value_serial, total_weight, serial_total_weight, nused, nused_serial;
    double seconds;
    
    // Local arrays
    int* scratch;
    int* scratch_cnt;
    
    // These constants have little effect on runtime
    int max_value  = 1000;
    int max_weight = 1000;
    
    // Shared memory arrays
    sintptr_nblk s_weight;
    sintptr_nblk s_value;
    sintptr_cblk s_old_row;
    sintptr_cblk s_new_row;
    sintptr_cblk s_old_cnt;
    sintptr_cblk s_new_cnt;
    sintptr_sngl s_progress;

    // Local pointers to shared memory
    int *l_weight;
    int *l_value;
    int *l_old_row;
    int *l_new_row;
    int *l_old_cnt;
    int *l_new_cnt;

    // Allocate shared memory arrays
    // Separate copies of weight and value are provided for each thread
    s_weight    = (sintptr_nblk) upc_all_alloc( THREADS, NITEMS * sizeof(int) );
    s_value     = (sintptr_nblk) upc_all_alloc( THREADS, NITEMS * sizeof(int) );
    s_old_row   = (sintptr_cblk) upc_all_alloc( THREADS, BLK_WIDTH * sizeof(int) );
    s_new_row   = (sintptr_cblk) upc_all_alloc( THREADS, BLK_WIDTH * sizeof(int) );
    s_old_cnt   = (sintptr_cblk) upc_all_alloc( THREADS, BLK_WIDTH * sizeof(int) );
    s_new_cnt   = (sintptr_cblk) upc_all_alloc( THREADS, BLK_WIDTH * sizeof(int) );
    s_progress  = (sintptr_sngl) upc_all_alloc( THREADS, sizeof(int) );

    if( !s_weight || !s_value || !s_old_row || !s_new_row || !s_old_cnt || !s_new_cnt || !s_progress)
    {
        fprintf( stderr, "Failed to allocate shared memory\n" );
        upc_global_exit( -1 );
    }

    // Allocate memory for local scratch arrays
    scratch     = malloc( BLK_WIDTH * sizeof(int) );
    scratch_cnt = malloc( BLK_WIDTH * sizeof(int) );
    
    if( !scratch || !scratch_cnt )
    {
        fprintf( stderr, "Failed to allocate local memory\n" );
        upc_global_exit( -1 );
    }
  
    // Initialize local pointers to shared memory
    // These pointers point to the processors own block in the shared memory arrays
    l_weight    = (int*)( s_weight + (MYTHREAD * NITEMS) );
    l_value     = (int*)( s_value + (MYTHREAD * NITEMS) );
    l_old_row   = (int*)( s_old_row + (MYTHREAD * BLK_WIDTH) );
    l_new_row   = (int*)( s_new_row + (MYTHREAD * BLK_WIDTH) );
    l_old_cnt   = (int*)( s_old_cnt + (MYTHREAD * BLK_WIDTH) );
    l_new_cnt   = (int*)( s_new_cnt + (MYTHREAD * BLK_WIDTH) );
    
    // Random initialization
    srand48( (unsigned int)time(NULL) );
    max_weight = min( max_weight, CAPACITY );//don't generate items that don't fit into bag
    
    // Only thread 0 must generate the random weights and values
    if( MYTHREAD == 0 )
    {
        for(i = 0; i < NITEMS; i++)
        {
            l_weight[i] = 1 + (lrand48()%max_weight);
            l_value[i]  = 1 + (lrand48()%max_value);
        }
    }

    upc_barrier;

    // Now copy the randomixed weights and values to the different processors
    if( MYTHREAD != 0 )
    {
        // Iterate through items
        for(i = 0; i < NITEMS; i++)
        {
            l_weight[i] = s_weight[i];
            l_value[i]  = s_value[i];
        }
    }
    
    upc_barrier;

    // time the solution
    seconds = read_timer();
    
    // Actually solve the knapsack problem
    // Must return best_value, and populate s_count array correctly
    best_value = solve_upc(NITEMS, CAPACITY, s_old_row, s_new_row, s_old_cnt, s_new_cnt, l_old_row, l_new_row,
                            l_old_cnt, l_new_cnt, scratch, scratch_cnt, l_weight, l_value, &nused,
                            &total_weight, s_progress);

    seconds = read_timer() - seconds;
    
    // Check the result against the serial code
    if( MYTHREAD == 0 )
    {
        // Print problems parameters
        printf( "%d items, capacity: %d, time: %g\n", NITEMS, CAPACITY, seconds );
        
        compute_serial(NITEMS, CAPACITY, l_weight, l_value, &nused_serial, &serial_total_weight, &best_value_serial);
        
        printf("Serial: %d items used, value %d, weight %d\n", nused, best_value_serial, serial_total_weight );
        
        // Print UPC solution
        printf( "%d items used, value %d, weight %d\n", nused, best_value, total_weight );
        
        // If different, print error message
        if( best_value != best_value_serial || nused != nused_serial || total_weight != serial_total_weight )
            printf( "WRONG SOLUTION\n" );
    

        // Free up allocated shared memory
        upc_free(s_weight);
        upc_free(s_value);
        upc_free(s_old_row);
        upc_free(s_new_row);
        upc_free(s_old_cnt);
        upc_free(s_new_cnt);
    }
    
    // Free up allocated local memory
    free(scratch);
    free(scratch_cnt);

    return 0;
}
