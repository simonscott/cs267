#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <upc.h>


//
// Constant Definitions
//
#define CAPACITY    999
#define NITEMS      5000
//#define CAPACITY    40
//#define NITEMS      8
#define BLK_WIDTH   (CAPACITY+1+THREADS-1) / THREADS

//
// Type Definitions
//    
typedef shared [BLK_WIDTH]  int* sintptr_cblk;
typedef shared [NITEMS]     int* sintptr_nblk;


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
int solve_serial( int nitems, int cap, int *w, int *v )
{
    int i, j, best, *allocated, *T, wj, vj;
    
    //alloc local resources
    T = allocated = malloc( nitems*(cap+1)*sizeof(int) );
    if( !allocated )
    {
        fprintf( stderr, "Failed to allocate memory" );
        upc_global_exit( -1 );
    }
    
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
    
    //free resources
    free( allocated );
    
    return best;
}


//
// Solves the knapsack problem using shared memory
//
int solve_upc(  int nitems, int cap, sintptr_cblk s_old_row, sintptr_cblk s_new_row, sintptr_cblk s_count, int* l_old_row, int* l_new_row, int* l_count,
                int* scratch, int* scratch_cnt, int* l_weight, int* l_value)
{
    // Local variables
    int i, j;
    int start_col, scratch_start_col, scratch_end_col;
    int *temp_intptr;
    sintptr_cblk temp_sintptr;

    // Intialize first row
    int val_item0 = l_value[0];

    upc_forall(i = 0; i < l_weight[0]; i++; &s_old_row[i])
        s_old_row[i] = 0;

    upc_forall(i = l_weight[0]; i < cap+1; i++; &s_old_row[i])
        s_old_row[i] = val_item0;

    upc_barrier;

    // Determine this thread's starting column
    start_col = MYTHREAD * BLK_WIDTH;

    // Iterate through remaining rows
    for(i = 1; i < nitems; i++)
    {
        // Copy the section of the above row, that we need and that belongs to another thread, to scratch
        // Note: before using scratch, took 15 seconds
        scratch_start_col = max(0, start_col - l_weight[i]);
        scratch_end_col = max(0, start_col - l_weight[i] + BLK_WIDTH);  // end is exclusive
        //upc_memget(scratch, &s_old_row[scratch_start_col], scratch_end_col - scratch_start_col);

        for(j = scratch_start_col; j < scratch_end_col; j++)
            scratch[j - scratch_start_col] = s_old_row[j];

        // Iterate through all columns belonging to this thread
        // Note: the last thread does some extra work on non-existent columns
        for(j = 0; j < BLK_WIDTH; j++)
        {
            // If this item is larger than the current capacity
            if(start_col + j < l_weight[i])
                l_new_row[j] = l_old_row[j];

            else
            {
                l_new_row[j] = max( l_old_row[j],
                                    scratch[start_col + j - l_weight[i] - scratch_start_col] + l_value[i] );
            }
        }

        // Sync all threads
        upc_barrier;

        // Swap new and old rows
        temp_intptr = l_old_row;
        l_old_row = l_new_row;
        l_new_row = temp_intptr;

        temp_sintptr = s_old_row;
        s_old_row = s_new_row;
        s_new_row = temp_sintptr;

        // Sync all threads (is this needed?)
        upc_barrier;
    }

    return s_old_row[cap]; //can we optimize this access?
}


//
//  benchmarking program
//
int main( int argc, char** argv )
{
    // Local variables
    int i, best_value, best_value_serial, total_weight, nused, total_value;
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
    sintptr_cblk s_count;

    // Local pointers to shared memory
    int *l_weight;
    int *l_value;
    int *l_old_row;
    int *l_new_row;
    int *l_count;

    // Allocate shared memory arrays
    // Separate copies of weight and value are provided for each thread
    s_weight    = (sintptr_nblk) upc_all_alloc( THREADS, NITEMS * sizeof(int) );
    s_value     = (sintptr_nblk) upc_all_alloc( THREADS, NITEMS * sizeof(int) );
    s_old_row   = (sintptr_cblk) upc_all_alloc( THREADS, BLK_WIDTH * sizeof(int) );
    s_new_row   = (sintptr_cblk) upc_all_alloc( THREADS, BLK_WIDTH * sizeof(int) );
    s_count     = (sintptr_cblk) upc_all_alloc( THREADS, BLK_WIDTH * sizeof(int) );

    if( !s_weight || !s_value || !s_old_row || !s_new_row || !s_count )
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
    l_count     = (int*)( s_count + (MYTHREAD * BLK_WIDTH) );
    
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
    
    // time the solution
    seconds = read_timer();
    
    // Actually solve the knapsack problem
    // Must return best_value, and populate s_count array correctly
    best_value = solve_upc(NITEMS, CAPACITY, s_old_row, s_new_row, s_count, l_old_row, l_new_row, l_count, scratch, scratch_cnt, l_weight, l_value);

    seconds = read_timer() - seconds;
    
    // Check the result against the serial code
    if( MYTHREAD == 0 )
    {
        // Debug output
        for(int i = 0; i < CAPACITY+1; i++)
        {
            printf("[%d] owned by THREAD %d\n", i, upc_threadof(&s_count[i]));
        }

        // Print problems parameters
        printf( "%d items, capacity: %d, time: %g\n", NITEMS, CAPACITY, seconds );
        
        best_value_serial = solve_serial( NITEMS, CAPACITY, l_weight, l_value );
        
        printf("Serial best value=%d. UPC best value=%d\n", best_value_serial, best_value);

        // Determine total weight and value, and num items used, by iterating through UPC output array
        total_weight = nused = total_value = 0;
        for( i = 0; i < NITEMS; i++ )
        {
            if( s_count[i] != 0 )
            {
                nused++;
                total_weight += s_weight[i];
                total_value += s_value[i];
            }
        }
        
        // Print UPC solution
        printf( "%d items used, value %d, weight %d\n", nused, total_value, total_weight );
        
        // If different, print error message
        if( best_value != best_value_serial || best_value != total_value || total_weight > CAPACITY )
            printf( "WRONG SOLUTION\n" );
    
        // Free up allocated shared memory
/*        upc_free(s_weight);
        upc_free(s_value);
        upc_free(s_old_row);
        upc_free(s_new_row);
        upc_free(s_count);
*/
    }
    
    // Free up allocated local memory
    free(scratch);
    free(scratch_cnt);

    return 0;
}
