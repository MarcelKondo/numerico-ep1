#include<stdio.h>
#include<math.h>
#include<vector>
#include<iostream>

#define EPS 1e-6


using namespace std;

// ----------Handler functions------------
bool isZero(double num){
    if(abs(num) < EPS) return true;
    else return false;
}
void printM(vector<vector<double>> &W, int n, int m){

    for(int i = 0 ; i < n ; i++){
        for( int j = 0 ; j < m ; j++){

            printf("% 1.5f", W[i][j]);
        }
        printf("\n");
    }
    printf("\n");

}

//---------------------------------------

void scCalc(double& s, double& c, double wjk, double wik){
    double tau;
    if( abs(wik) > abs(wjk) ){
        tau = (-1.0)*wjk/wik;
        c = 1/sqrt(1 + tau*tau);
        s = c*tau;
    }
    else{
        tau = (-1.0)*wik/wjk;
        s = 1/sqrt(1 + tau*tau);
        c = s*tau;
    }
    
    return;
}

void RotGivens(vector<vector<double>> &W, int n, int m, int i, int j, double s, double c){
    double aux;
    for(int r = 0 ; r < m ; r++ ){

        aux = c*W[i][r] - s*W[j][r];
        W[j][r] = s*W[i][r] + c*W[j][r];
        W[i][r] = aux;

    }
}

void SolveSys(vector<vector<double>> &W, int n, int m, vector<vector<double>> &b, vector<vector<double>> &x){
    int i, j, k;
    double s, c;
    for( k = 0 ; k < m ; k++){
        for( j = n-1 ; j >= k+1 ; j--){
            i = j-1;
            if( !isZero(W[j][k]) ){
                scCalc( s, c, W[j][k], W[i][k]);
                RotGivens( W, n, m, i, j, s, c);
                RotGivens( b, m, b[0].size(), i, j, s, c);
            }
        }
    }

    for(int p = 0 ; p < b[0].size() ; p++){ //solving systems in parallel
        for( k = m-1 ; k >= 0 ; k--){
            double aux = 0.0;
            
            for( j = k + 1 ;  j < m ; j++) aux += W[k][j] * x[j][p];
            
            x[k][p] = (b[k][p] - aux)/W[k][k];
           
        }
       
    }

}


int main(){
    int n, m;
    vector<vector<double>> W;
    vector<vector<double>> b;
    vector<vector<double>> x;
    //Task1
    {
        //item a
        printf("START 1-A\n");
        n = 10;
        m = 10;
        W.assign(n, vector<double>(m, 0.0));
        b.assign(n, vector<double>(1, 1.0));
        x.assign(m, vector<double>(1, 0.0));
        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < m ; j++){
                if(i == j) W[i][j] = 2.0;
                else if(abs(i - j) == 1) W[i][j] = 1.0;
                else W[i][j] = 0.0;
            }
        }
        printM(W,n,m);
        SolveSys(W, n, m, b, x);
        printM(x, m, 1);

        //end a
        //item b
        printf("START 1-B\n");
        n = 20;
        m = 17;
        W.assign(n, vector<double>(m, 0.0));
        b.assign(n, vector<double>(1, 1.0));
        x.assign(m, vector<double>(1, 0.0));
        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < m ; j++){
                if(abs(i-j) <= 4) W[i][j] = 1.0/(i+1 + j+1 -1);
                else W[i][j] = 0.0;
            }
            b[i][0] = i+1;
        }
        printM(W,n,m);
        SolveSys(W, n, m, b, x);
        printM(x, m, 1);
        //end b

         //item c
        printf("START 1-C\n");
        n = 10;
        m = 10;
        W.assign(n, vector<double>(m, 0.0));
        b.assign(n, vector<double>(3, 1.0));
        x.assign(m, vector<double>(3, 0.0));
        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < m ; j++){
                if(i == j) W[i][j] = 2.0;
                else if(abs(i - j) == 1) W[i][j] = 1.0;
                else W[i][j] = 0.0;
            }
            b[i][1] = i+1;
            b[i][2] = 2*(i+1) - 1;
        }
        printM(W,n,m);
        SolveSys(W, n, m, b, x);
        printM(x, m, 3);

        //end c
        //item d
        printf("START 1-D\n");
        n = 20;
        m = 17;
        W.assign(n, vector<double>(m, 0.0));
        b.assign(n, vector<double>(3, 1.0));
        x.assign(m, vector<double>(3, 0.0));
        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < m ; j++){
                if(abs(i-j) <= 4) W[i][j] = 1.0/(i+1 + j+1 -1);
                else W[i][j] = 0.0;
            }
            b[i][1] = i+1;
            b[i][2] = 2*(i+1) - 1;
        }
        printM(W,n,m);
        SolveSys(W, n, m, b, x);
        
        printM(x, m, 3);
        //end d
        
    }

    return 0;
}
