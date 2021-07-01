//
//  main.cpp
//  nCr
//
//  Created by Hiren Thummar on 8/8/19.
//  Copyright © 2019 Hiren Thummar. All rights reserved.
//
// C++ program to print all combination of size r in an array of size n
//#include <bits/stdc++.h>
//

#include <iostream>
#include <bitset>

using namespace std;

void combinationUtil(int arr[], int data[],
                     int start, int end,
                     int index, int r);

// The main function that prints all combinations of size r in arr[] of size n.
// This function mainly uses combinationUtil()
void printCombination(int arr[], int n, int r)
{
    // A temporary array to store all combination one by one
    int data[r];
    
    // Print all combination using temprary array 'data[]'
    combinationUtil(arr, data, 0, n-1, 0, r);
}

/* arr[] ---> Input Array
 data[] ---> Temporary array to
 store current combination
 start & end ---> Staring and
 Ending indexes in arr[]
 index ---> Current index in data[]
 r ---> Size of a combination to be printed */
void combinationUtil(int arr[], int data[],
                     int start, int end,
                     int index, int r)
{
    // Current combination is ready
    // to be printed, print it
    if (index == r) {
        for (int j = 0; j < r; j++) {
            cout << data[j] << " ";
        }
        cout << endl;
    }
    
    // replace index with all possible
    // elements. The condition "end-i+1 >= r-index"
    // makes sure that including one element
    // at index will make a combination with
    // remaining elements at remaining positions
    for (int i = start; i <= end &&
         end - i + 1 >= r - index; i++)
    {
        data[index] = arr[i];
        combinationUtil(arr, data, i+1,
                        end, index+1, r);
    }
}

// Driver code
int main() {
    int arr[] = {1, 2, 3, 4, 5, 6, 7, 8};
//    int r = 3;
    int n = sizeof(arr)/sizeof(arr[0]);
    for (int i = 0; i < 9; ++i){
        printCombination(arr, n, i);
    }
}
