/*
    Classic McEliece Parameters: (m, t)
    1. (12, 64)
    2. (13, 96)
    3. (13, 128)
    4. (13, 119)
    5. (13, 128)
*/

/**********************
 * 사용해야하는 파라미터 
 ***********************/

#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include "fq_nmod.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "fq_nmod_poly.h"
#include "nmod_poly.h"
#include "fmpz.h"

#define UPBOUND 10000  // 2^m -1 개의 테이블 생성 
#define FALSE 0    // 두 값을 비교하는 과정에서 false인 경우 0을 출력함.

    int
    main()
{
    double t_precomp=0, t_modexp=0;
    
    int n_iter=10;

    for(int cnt=0;cnt<n_iter;cnt++){

    fmpz_t q, twom, twom_1, twomt_1, two;
    slong  m, z, t;
    ulong  u1, u2, u3;
    clock_t c0, c1;
    fq_nmod_ctx_t ctx;
    fq_nmod_t qp[UPBOUND], qi[UPBOUND], coef, table[UPBOUND];
    fq_nmod_poly_t Qx, g, Sq, x, r[UPBOUND], result, qod, qev;
    nmod_poly_t np1, np2;

    FLINT_TEST_INIT(state);    // 랜덤 생성을 위한 시드 느낌.
    
    //============다항식 구조체 생성========================================
    fq_nmod_poly_init(Qx, ctx);  // Q(x)
    fq_nmod_poly_init(g, ctx);   // g(x)
    fq_nmod_poly_init(Sq, ctx);  // Sq(x)
    fq_nmod_poly_init(x, ctx);   // x(x) = x
    //============사용자가 설정하는 변수===============================    
    fmpz_init_set_ui(q, 2);       // Fqm에서 q
    m = 12;                       // Fqm에서 m
    t = 41;                       // Fqm(X)의 차수 설정. g(x)는 t, Q(x)는 t-1차 다항식.
    //============================================================    
    fq_nmod_ctx_init_conway(ctx, q, m, "t");                   // ctx에 Fqm저장. 
    //printf("\n======ctx info======\n");    fq_nmod_ctx_print(ctx);   printf("====================\n"); // 저장된 ctx정보 확인
    fq_nmod_poly_gen(x, ctx);                                  // 1차다항식 x(x)=x설정. 
    // 다항식 Qx,g는 랜덤하게 생성.    
    fq_nmod_poly_randtest(Qx, state, t, ctx);                  // t개의 항을 가진 Qx. 즉 t-1차가 됨.
    fq_nmod_poly_randtest_irreducible(g, state, t+1, ctx);     // t+1개의 항을 가진 g. 즉 t차가 됨. 단 g는 기약다항식으로 설정.
    //==================== pre-computation ===========================================================
    //=========================== Ri(X) 사전계산 LINE 2 ==================================
    fmpz_init(two);               // 2 저장
    fmpz_init(twomt_1);           // 2^(mt-1) 저장
    fmpz_set_ui(two, 2);
    fmpz_pow_ui(twomt_1, two, (m*t-1));
    fq_nmod_poly_powmod_fmpz_binexp(Sq, x, twomt_1, g, ctx);    // Sq(x) = x^(2^(mt-1)) mod g(x)
    //printf("Sq(x) = ");     fq_nmod_poly_print_pretty(Sq, "X", ctx);    printf("\n");   // Sq 확인   
    for(ulong i = 0; i <= (t/2 - 1); i++)
    {
        fq_nmod_poly_init(r[i], ctx);                           // r[i] 구조체 생성.
        fq_nmod_poly_powmod_ui_binexp(r[i], x, i, g, ctx);      // r[i] = x^i mod g(x)
        fq_nmod_poly_mulmod(r[i], r[i], Sq, g, ctx);            // r[i] = r[i]*Sq mod g(x)
        //printf("ri(x) = ");     fq_nmod_poly_print_pretty(r[i],"X",ctx);    printf("\n");   //ri 확인 
    }
    //======== 사전계산테이블(Fqm의 모든 원소에 대해 2^(m-1)제곱을 미리 해놓은 테이블) LINE 3 =========
    fmpz_init(twom_1);                                      // twom_1 = 2^(m-1). (fmpz)
    fmpz_init(twom);                                        // twom = 2^m        (fmpz)
    fmpz_set_ui(twom_1, pow(2,m-1));                        
    fmpz_set_ui(twom, pow(2,m));
    u3 = fmpz_get_ui(twom);                                 // u3 = 2^m (ulong)
    //printf("2^(m-1) = ");  fmpz_print(twom_1);  printf("\n");
    
    fq_nmod_init(coef,ctx);                                  
    for(ulong i = 0; i < u3; i++)
    {    
        nmod_poly_init(np1, m);      // 차수는 m차 미만. 
        fq_nmod_init(table[i], ctx); // fqm상의 모든 원소에 대해 2^twom_1 을 table에 저장.
        u1 = i;
        for(ulong j = 0; j < m; j++) // 해당 For문은 주어진 j라는 정수를 유한체 2^m상의 원소로 바꾸는 과정을 거친다. 
        {
            u2 = u1 % 2;
            u1 = u1 / 2;
            nmod_poly_set_coeff_ui(np1, j, u2); // 2씩 나누면서 fqm의 원소로 표현.
        }
        // ========테이블 생성 확인=========
        // printf("%lu  ", i);                                // i번째 테이블
        // nmod_poly_print_pretty(np1, "x");  printf("  ");   // nmod 형태 출력 
        fq_nmod_set_nmod_poly(coef, np1, ctx);                // nmod를 fqnmod 형태로
        fq_nmod_pow(table[i], coef, twom_1 ,ctx);             // table = coef^twom_1 
        // fq_nmod_print_pretty(coef, ctx); printf(" ");      // fqnmod 형태 출력 
        // fq_nmod_print_pretty(table[i], ctx); printf("\n"); // table 출력
        nmod_poly_clear(np1);
    }
    //printf("\n=========\n");
    //================================================================================================

    //printf("x    = ");     fq_nmod_poly_print_pretty(x, "X", ctx);     printf("\n");  // x확인
    //printf("g(x) = ");     fq_nmod_poly_print_pretty(g, "X", ctx);     printf("\n");  // g확인
    //printf("Q(x) = ");     fq_nmod_poly_print_pretty(Qx, "X", ctx);    printf("\n");  // Q확인
    
    nmod_poly_init(np2, t);  
    for(slong i=0; i < t; i++)
    {
        fq_nmod_init(qi[i], ctx);                          // fqm원소 qi[i]
        fq_nmod_init(qp[i], ctx); 
    }
    fq_nmod_poly_init(result, ctx); 
    fq_nmod_poly_init(qev, ctx); fq_nmod_poly_init(qod, ctx);  
    
    int pow2[20]={0,};
    for(int j=0;j<m;j++)
    {
        pow2[j]=pow(2, j);
    }

    //============= LINE1 ==========================  
    c0 = clock();                                            // time start   
    // LINE 1는 LINE 2를 위함이나 LINE 2는 사전계산으로 진행했으므로 진행할 필요없음.  
    // printf("Sq(x) = ");     fq_nmod_poly_print_pretty(Sq, "X", ctx);    printf("\n");  //Sq 확인

    //============= LINE2 ==========================     
    // LINE 2는 사전계산을 통해 진행. 
    
    //============= LINE3 ==========================     
    for(slong i=0; i < t; i++)
    {
        //nmod_poly_init(np2, t);                            // t-1차 이하 -> 초기화 
        nmod_poly_zero(np2);                               // t-1차 이하 -> 초기화 
        //fq_nmod_init(qi[i], ctx);                          // fqm원소 qi[i]
        fq_nmod_poly_get_coeff(qi[i], Qx, i, ctx);         // 입력받은 Qx의 계수를 추출. 현재는 fqm원소 형태 
        fq_nmod_get_nmod_poly(np2, qi[i], ctx);            // 해당 원소를 x이용 방정식 형태로 변환
        //printf("coef = "); nmod_poly_print_pretty(np2, "x");
        
        z = 0;                                             // 위의 Fqm원소를 숫자로 바꾸어 저장.  
        for(slong j = 0; j < m; j++)                       // 해당 for문은 Fqm원소를 숫자로 바꾸어 저장.
        { 
            z += nmod_poly_get_coeff_ui(np2, j) * pow2[j];                   // 숫자형태로 변환
        }
        //==============결과 확인 ===============================
        //printf(" z=%ld ",z);                                                   // 위의 x형태의 정수형태
        //printf("table= "); fq_nmod_print_pretty(table[z], ctx); printf("\n");  // 해당테이블에 저장된값
        //fq_nmod_init(qp[i], ctx); 
        fq_nmod_set(qp[i], table[z], ctx);                                     // 사전계산한 값을 qp에 저장
        //printf("qp= "); fq_nmod_print_pretty(qp[i], ctx); 

        //nmod_poly_print_pretty(np2, "x"); printf("\n");
        //nmod_poly_clear(np2);
    }
    
    //============= LINE4, 5, 6 ==========================     
    //fq_nmod_poly_init(result, ctx);   // 최종 결과 저장 변수
    fq_nmod_poly_zero(result, ctx);   // result를 0으로 초기화

    for(int i = 0; i < t; i++)
    {
        //printf("qp= "); fq_nmod_print_pretty(qp[i], ctx); printf("  ");            // 값확인
        if(i % 2 == 0)  // index가 짝수인 경우 LINE 4의 연산 따름
        {
            //fq_nmod_poly_init(qev, ctx);   
            fq_nmod_poly_zero(qev,ctx);                                          // qp의 index가 짝수인 경우 중간값 저장에 사용
            fq_nmod_poly_set_coeff(qev, i/2, qp[i], ctx);                            // qev를 qp[i]x^(i/2)로 설정.
            fq_nmod_poly_add(result, result, qev, ctx);                              // 해당 식을 최종결과에 더함
            //printf("q'= "); fq_nmod_print_pretty(qp[i], ctx);                        // 결과확인
            //printf("qev= "); fq_nmod_poly_print_pretty(qev, "X", ctx); printf("  "); // 결과확인
            //fq_nmod_poly_clear(qev, ctx);                                            // 변수 초기화
        }
        else // index가 홀수인 경우 LINE 5의 연산 따름
        {    
            //fq_nmod_poly_init(qod, ctx);                                                  // qp의 index가 홀수인 경우 중간값 저장에 사용
            fq_nmod_poly_zero(qod,ctx); 
            fq_nmod_poly_scalar_mul_fq_nmod(qod, r[i/2], qp[i], ctx);                     // qod를 r[i/2]에 qp[i]배 한 다항식으로 설정 
            fq_nmod_poly_add(result, result, qod, ctx);                                   // 해당 식을 최종결과에 더함    
            //printf("q'= "); fq_nmod_print_pretty(qp[i], ctx);                             // 결과확인
            //printf("r[i]= "); fq_nmod_poly_print_pretty(r[i/2], "X", ctx); printf("  ");  // 결과확인
            //printf("qr[i]="); fq_nmod_poly_print_pretty(qod, "X", ctx);    printf("  ");  // 결과확인
            //fq_nmod_poly_clear(qod, ctx);                                            
        }
        //printf("return(x) = "); fq_nmod_poly_print_pretty(result, "X", ctx);  printf("\n");   // ri 
    }
    c1 = clock();                           // time end
    t_precomp += (double) (c1 - c0) / CLOCKS_PER_SEC;

    nmod_poly_clear(np2);
    fq_nmod_poly_clear(qev, ctx);     // 변수 초기화
    fq_nmod_poly_clear(qod, ctx);                                            


    //printf("\n\n");
    //printf("g(x)      = ");    fq_nmod_poly_print_pretty(g, "X", ctx);      printf("\n\n");  // g확인
    //printf("Q(x)      = ");    fq_nmod_poly_print_pretty(Qx, "X", ctx);     printf("\n\n");  // Q확인
    //printf("return(x) = ");    fq_nmod_poly_print_pretty(result, "X", ctx); printf("\n\n");  // sqrt(Qx) 출력
    
    fq_nmod_poly_t calre;                               // 계산값이 실제 값과 같은지 비교하기 위함.  최종결과인 result의 제곱이 Q(x)가 나오는지 확인
    
    fq_nmod_poly_init(calre, ctx);
    fq_nmod_poly_mulmod(calre, result, result, g, ctx); //result의 제곱을 mod g(x)에 대해 계산하여 calre에 설정.
    
    //printf("cal       = ");    fq_nmod_poly_print_pretty(calre, "X", ctx);  printf("\n");  //ri

    int a,b;
    
    a = fq_nmod_poly_equal(Qx,calre,ctx);               // Qx == calre 인 경우 1출력, 아니면 0출력

    //printf("\n\n");
    //printf("m = %ld, t = %ld \n", m, t);
    //printf("%d",a);
    if(a==FALSE)
    {
        printf("\n FALSE \n");
        break;
    }
    else
        printf("\n TRUE Q(x) == return^2 \n");
    
    //printf("time algori : %fs\n", t_precomp);     // 알고리즘 계산시간

    //================== 지수승 연산 ================================
    // fqm[X]/g(X)의 원소 f(X)에 대해 f(X)^(2^(mt-1)) mod g(X)는 sqrt(f(X))가 됨을 이용.
    fq_nmod_poly_t expone;          
    fq_nmod_poly_init(expone,ctx);


    c0 = clock();                       // time start
    fq_nmod_poly_powmod_fmpz_binexp(expone, Qx, twomt_1, g, ctx);     // Qx^(2^(mt-1)) mod g(X) 
    //printf("Sq(x) = ");     fq_nmod_poly_print_pretty(h,"X",ctx);    printf("\n");  //Sq 확인
    c1 = clock();                       // time end
    
    
    t_modexp += (double) (c1 - c0) / CLOCKS_PER_SEC;

    b = fq_nmod_poly_equal(result,expone,ctx);          // 위의 알고리듬과 같은 결과가 나오는지 확인. 
    
    if(b==FALSE) //결과 확인용
    {
        printf("\n FALSE \n");
        break;
    }
    else
        printf(" TRUE sqrt(Q(x)) == return \n");   

    //printf("time expone : %fs\n", t_modexp); // 지수승 계산시간
    //printf("\n\n");
    
    //====================== 메모리 해제 ============================
    fmpz_clear(q);  fmpz_clear(twom);  fmpz_clear(twom_1);  fmpz_clear(twomt_1);  fmpz_clear(two);  

    for (int i = 0; i < t; i++)
    {
        fq_nmod_clear(qp[i],ctx);
        fq_nmod_clear(qi[i],ctx);
    }   
    for (int i = 0; i < u3; i++)
        fq_nmod_clear(table[i],ctx);
    for (int i = 0; i < t/2 ; i++)
        fq_nmod_poly_clear(r[i],ctx);
    
    fq_nmod_clear(coef,ctx);

    fq_nmod_poly_clear(Qx,ctx);
    fq_nmod_poly_clear(g,ctx);
    fq_nmod_poly_clear(Sq,ctx);
    fq_nmod_poly_clear(x,ctx);
    fq_nmod_poly_clear(result,ctx);
    fq_nmod_poly_clear(calre,ctx);
    fq_nmod_poly_clear(expone,ctx);

    fq_nmod_ctx_clear(ctx);
    printf("%d\n",cnt);
    
    }

    printf("\n");
    printf("time algori : %fs\n", t_precomp/n_iter);     // 알고리즘 계산시간
    printf("time expone : %fs\n", t_modexp/n_iter);     // 지수승 계산시간


    return EXIT_SUCCESS;
}