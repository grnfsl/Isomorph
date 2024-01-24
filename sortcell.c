/* 
 * This file contains an efficient sorting algorithm for sorting vertices
 * according to their degree and also check if cell contains distinct degrees
 * which is used to check if partition become discrete
 */

#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#ifdef _OPENMP
#include <omp.h>
#endif

const unsigned long cutoff_seq_par = 160000L; //cutoff for sequential and parallel quick sort
const unsigned long cutoff_ins_quk = 32L; //cutoff for insertion and quick sort

int verify_sort(int *arr, int *degree, int n)
{
    for(int i = 1; i < n; ++i)
        if(degree[arr[i]] < degree[arr[i-1]]) return 0;
    return 1;
}

static inline void swap(int *a, int *b)
{
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

static inline unsigned long min_ind(int *arr, int *degree, unsigned long p, unsigned long r) //find index of minimum element
{
    int min = arr[p];
    unsigned long ind = p;
    ++p;
    while(p < r+1){
        if(degree[arr[p]] < degree[min]){
            min = arr[p];
            ind = p;
        }
        ++p;
    }
    return ind;
}

void insertion_sort(int *arr, int *degree, unsigned long p, unsigned long r)
//constraints: p<=r
{
    int key;
    unsigned long i, j, ind;

    ind = min_ind(arr, degree, p, r);
    swap(&arr[p], &arr[ind]);

    for(i = p+2; i < r+1; ++i){
        key = arr[i];
        for (j = i; degree[key] < degree[arr[j-1]]; --j) {
            arr[j] = arr[j-1];
        }
        arr[j] = key;
    }
}

int median3(int *a, int *degree, unsigned long p, unsigned long r)
//arr must contain at least 3 elements and p<r and r>0
{
    int center;

    center = (p+r)/2;

    if(degree[a[p]] > degree[a[center]])
        swap(&a[p], &a[center]);
    if(degree[a[center]] > degree[a[r]]){
        swap(&a[center], &a[r]);
        if(degree[a[p]] > degree[a[center]])
            swap(&a[p], &a[center]);
    }
    swap(&a[center], &a[r-1]);
    return a[r-1];
}

void sq_sort(int *arr, int *degree, unsigned long p, unsigned long r)
{
    unsigned long i, j;
    int pivot;

    if((p+cutoff_ins_quk) > r){
        insertion_sort(arr, degree, p, r);
    }
    else{
        pivot = median3(arr, degree, p, r);
        i = p;
        j = r-1;
        while (1) {
            while (degree[arr[++i]] < degree[pivot]);
            while (degree[arr[--j]] > degree[pivot]);
            if(i < j)
                swap(&arr[i], &arr[j]);
            else break;
        }
        swap(&arr[i], &arr[r-1]);
        sq_sort(arr, degree, p, i-1);
        sq_sort(arr, degree, i+1, r);
    }
}

void pq_sort(int *arr, int *degree, unsigned long p, unsigned long r)
{
    unsigned long i, j;
    int pivot;

    if((p+cutoff_ins_quk) > r){
        insertion_sort(arr, degree, p, r);
    }
    else{
        pivot = median3(arr, degree, p, r);
        i = p;
        j = r-1;
        while (1) {
            while (degree[arr[++i]] < degree[pivot]);
            while (degree[arr[--j]] > degree[pivot]);
            if(i < j){
                if(degree[arr[i]] != degree[arr[j]])
                    swap(&arr[i], &arr[j]);
            }
            else break;
        }
        swap(&arr[i], &arr[r-1]);

        #pragma omp task
        pq_sort(arr, degree, p, i-1);
        #pragma omp task
        pq_sort(arr, degree, i+1, r);
    }
}

bool qsortcell(int *arr, int *degree, unsigned long n)
{
    if(n < cutoff_seq_par){
        sq_sort(arr, degree, 0, n-1);
    }
    else{
        #ifdef _OPENMP
        #pragma omp parallel num_threads(8)
        #pragma omp single
        pq_sort(arr, degree, 0, n-1);
        #else
        sq_sort(arr, degree, 0, n-1);
        #endif
    }

    unsigned long j;
    for(j = 1; j < n; ++j)
        if(degree[arr[j]] == degree[arr[j-1]]) return false;
    return true;
}





































