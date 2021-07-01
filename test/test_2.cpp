//
//  test_2.cpp
//  hwlqmc_test
//
//  Created by Hiren Thummar on 2/26/21.
//  Copyright © 2021 Hiren Thummar. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;
int main(){
    srand((unsigned int) time(NULL));
    
    for(int i=1; i<=28; ++i){
        cout << i << ": " <<  rand() % 28 << endl;
        //    cout << rand() << endl;
    }
    return 0;
    
}
