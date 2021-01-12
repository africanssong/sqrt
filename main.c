/*
    Classic McEliece Parameters: (m, t)
    1. (12, 64)
    2. (13, 96)
    3. (13, 128)
    4. (13, 119)
    5. (13, 128)
*/

/**********************
 * ����ؾ��ϴ� �Ķ���� 
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

#define UPBOUND 10000  // 2^m -1 ���� ���̺� ���� 
#define FALSE 0    // �� ���� ���ϴ� �������� false�� ��� 0�� �����.

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

    FLINT_TEST_INIT(state);    // ���� ������ ���� �õ� ����.
    
    //============���׽� ����ü ����========================================
    fq_nmod_poly_init(Qx, ctx);  // Q(x)
    fq_nmod_poly_init(g, ctx);   // g(x)
    fq_nmod_poly_init(Sq, ctx);  // Sq(x)
    fq_nmod_poly_init(x, ctx);   // x(x) = x
    //============����ڰ� �����ϴ� ����===============================    
    fmpz_init_set_ui(q, 2);       // Fqm���� q
    m = 12;                       // Fqm���� m
    t = 41;                       // Fqm(X)�� ���� ����. g(x)�� t, Q(x)�� t-1�� ���׽�.
    //============================================================    
    fq_nmod_ctx_init_conway(ctx, q, m, "t");                   // ctx�� Fqm����. 
    //printf("\n======ctx info======\n");    fq_nmod_ctx_print(ctx);   printf("====================\n"); // ����� ctx���� Ȯ��
    fq_nmod_poly_gen(x, ctx);                                  // 1�����׽� x(x)=x����. 
    // ���׽� Qx,g�� �����ϰ� ����.    
    fq_nmod_poly_randtest(Qx, state, t, ctx);                  // t���� ���� ���� Qx. �� t-1���� ��.
    fq_nmod_poly_randtest_irreducible(g, state, t+1, ctx);     // t+1���� ���� ���� g. �� t���� ��. �� g�� �����׽����� ����.
    //==================== pre-computation ===========================================================
    //=========================== Ri(X) ������� LINE 2 ==================================
    fmpz_init(two);               // 2 ����
    fmpz_init(twomt_1);           // 2^(mt-1) ����
    fmpz_set_ui(two, 2);
    fmpz_pow_ui(twomt_1, two, (m*t-1));
    fq_nmod_poly_powmod_fmpz_binexp(Sq, x, twomt_1, g, ctx);    // Sq(x) = x^(2^(mt-1)) mod g(x)
    //printf("Sq(x) = ");     fq_nmod_poly_print_pretty(Sq, "X", ctx);    printf("\n");   // Sq Ȯ��   
    for(ulong i = 0; i <= (t/2 - 1); i++)
    {
        fq_nmod_poly_init(r[i], ctx);                           // r[i] ����ü ����.
        fq_nmod_poly_powmod_ui_binexp(r[i], x, i, g, ctx);      // r[i] = x^i mod g(x)
        fq_nmod_poly_mulmod(r[i], r[i], Sq, g, ctx);            // r[i] = r[i]*Sq mod g(x)
        //printf("ri(x) = ");     fq_nmod_poly_print_pretty(r[i],"X",ctx);    printf("\n");   //ri Ȯ�� 
    }
    //======== ����������̺�(Fqm�� ��� ���ҿ� ���� 2^(m-1)������ �̸� �س��� ���̺�) LINE 3 =========
    fmpz_init(twom_1);                                      // twom_1 = 2^(m-1). (fmpz)
    fmpz_init(twom);                                        // twom = 2^m        (fmpz)
    fmpz_set_ui(twom_1, pow(2,m-1));                        
    fmpz_set_ui(twom, pow(2,m));
    u3 = fmpz_get_ui(twom);                                 // u3 = 2^m (ulong)
    //printf("2^(m-1) = ");  fmpz_print(twom_1);  printf("\n");
    
    fq_nmod_init(coef,ctx);                                  
    for(ulong i = 0; i < u3; i++)
    {    
        nmod_poly_init(np1, m);      // ������ m�� �̸�. 
        fq_nmod_init(table[i], ctx); // fqm���� ��� ���ҿ� ���� 2^twom_1 �� table�� ����.
        u1 = i;
        for(ulong j = 0; j < m; j++) // �ش� For���� �־��� j��� ������ ����ü 2^m���� ���ҷ� �ٲٴ� ������ ��ģ��. 
        {
            u2 = u1 % 2;
            u1 = u1 / 2;
            nmod_poly_set_coeff_ui(np1, j, u2); // 2�� �����鼭 fqm�� ���ҷ� ǥ��.
        }
        // ========���̺� ���� Ȯ��=========
        // printf("%lu  ", i);                                // i��° ���̺�
        // nmod_poly_print_pretty(np1, "x");  printf("  ");   // nmod ���� ��� 
        fq_nmod_set_nmod_poly(coef, np1, ctx);                // nmod�� fqnmod ���·�
        fq_nmod_pow(table[i], coef, twom_1 ,ctx);             // table = coef^twom_1 
        // fq_nmod_print_pretty(coef, ctx); printf(" ");      // fqnmod ���� ��� 
        // fq_nmod_print_pretty(table[i], ctx); printf("\n"); // table ���
        nmod_poly_clear(np1);
    }
    //printf("\n=========\n");
    //================================================================================================

    //printf("x    = ");     fq_nmod_poly_print_pretty(x, "X", ctx);     printf("\n");  // xȮ��
    //printf("g(x) = ");     fq_nmod_poly_print_pretty(g, "X", ctx);     printf("\n");  // gȮ��
    //printf("Q(x) = ");     fq_nmod_poly_print_pretty(Qx, "X", ctx);    printf("\n");  // QȮ��
    
    nmod_poly_init(np2, t);  
    for(slong i=0; i < t; i++)
    {
        fq_nmod_init(qi[i], ctx);                          // fqm���� qi[i]
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
    // LINE 1�� LINE 2�� �����̳� LINE 2�� ����������� ���������Ƿ� ������ �ʿ����.  
    // printf("Sq(x) = ");     fq_nmod_poly_print_pretty(Sq, "X", ctx);    printf("\n");  //Sq Ȯ��

    //============= LINE2 ==========================     
    // LINE 2�� ��������� ���� ����. 
    
    //============= LINE3 ==========================     
    for(slong i=0; i < t; i++)
    {
        //nmod_poly_init(np2, t);                            // t-1�� ���� -> �ʱ�ȭ 
        nmod_poly_zero(np2);                               // t-1�� ���� -> �ʱ�ȭ 
        //fq_nmod_init(qi[i], ctx);                          // fqm���� qi[i]
        fq_nmod_poly_get_coeff(qi[i], Qx, i, ctx);         // �Է¹��� Qx�� ����� ����. ����� fqm���� ���� 
        fq_nmod_get_nmod_poly(np2, qi[i], ctx);            // �ش� ���Ҹ� x�̿� ������ ���·� ��ȯ
        //printf("coef = "); nmod_poly_print_pretty(np2, "x");
        
        z = 0;                                             // ���� Fqm���Ҹ� ���ڷ� �ٲپ� ����.  
        for(slong j = 0; j < m; j++)                       // �ش� for���� Fqm���Ҹ� ���ڷ� �ٲپ� ����.
        { 
            z += nmod_poly_get_coeff_ui(np2, j) * pow2[j];                   // �������·� ��ȯ
        }
        //==============��� Ȯ�� ===============================
        //printf(" z=%ld ",z);                                                   // ���� x������ ��������
        //printf("table= "); fq_nmod_print_pretty(table[z], ctx); printf("\n");  // �ش����̺� ����Ȱ�
        //fq_nmod_init(qp[i], ctx); 
        fq_nmod_set(qp[i], table[z], ctx);                                     // ��������� ���� qp�� ����
        //printf("qp= "); fq_nmod_print_pretty(qp[i], ctx); 

        //nmod_poly_print_pretty(np2, "x"); printf("\n");
        //nmod_poly_clear(np2);
    }
    
    //============= LINE4, 5, 6 ==========================     
    //fq_nmod_poly_init(result, ctx);   // ���� ��� ���� ����
    fq_nmod_poly_zero(result, ctx);   // result�� 0���� �ʱ�ȭ

    for(int i = 0; i < t; i++)
    {
        //printf("qp= "); fq_nmod_print_pretty(qp[i], ctx); printf("  ");            // ��Ȯ��
        if(i % 2 == 0)  // index�� ¦���� ��� LINE 4�� ���� ����
        {
            //fq_nmod_poly_init(qev, ctx);   
            fq_nmod_poly_zero(qev,ctx);                                          // qp�� index�� ¦���� ��� �߰��� ���忡 ���
            fq_nmod_poly_set_coeff(qev, i/2, qp[i], ctx);                            // qev�� qp[i]x^(i/2)�� ����.
            fq_nmod_poly_add(result, result, qev, ctx);                              // �ش� ���� ��������� ����
            //printf("q'= "); fq_nmod_print_pretty(qp[i], ctx);                        // ���Ȯ��
            //printf("qev= "); fq_nmod_poly_print_pretty(qev, "X", ctx); printf("  "); // ���Ȯ��
            //fq_nmod_poly_clear(qev, ctx);                                            // ���� �ʱ�ȭ
        }
        else // index�� Ȧ���� ��� LINE 5�� ���� ����
        {    
            //fq_nmod_poly_init(qod, ctx);                                                  // qp�� index�� Ȧ���� ��� �߰��� ���忡 ���
            fq_nmod_poly_zero(qod,ctx); 
            fq_nmod_poly_scalar_mul_fq_nmod(qod, r[i/2], qp[i], ctx);                     // qod�� r[i/2]�� qp[i]�� �� ���׽����� ���� 
            fq_nmod_poly_add(result, result, qod, ctx);                                   // �ش� ���� ��������� ����    
            //printf("q'= "); fq_nmod_print_pretty(qp[i], ctx);                             // ���Ȯ��
            //printf("r[i]= "); fq_nmod_poly_print_pretty(r[i/2], "X", ctx); printf("  ");  // ���Ȯ��
            //printf("qr[i]="); fq_nmod_poly_print_pretty(qod, "X", ctx);    printf("  ");  // ���Ȯ��
            //fq_nmod_poly_clear(qod, ctx);                                            
        }
        //printf("return(x) = "); fq_nmod_poly_print_pretty(result, "X", ctx);  printf("\n");   // ri 
    }
    c1 = clock();                           // time end
    t_precomp += (double) (c1 - c0) / CLOCKS_PER_SEC;

    nmod_poly_clear(np2);
    fq_nmod_poly_clear(qev, ctx);     // ���� �ʱ�ȭ
    fq_nmod_poly_clear(qod, ctx);                                            


    //printf("\n\n");
    //printf("g(x)      = ");    fq_nmod_poly_print_pretty(g, "X", ctx);      printf("\n\n");  // gȮ��
    //printf("Q(x)      = ");    fq_nmod_poly_print_pretty(Qx, "X", ctx);     printf("\n\n");  // QȮ��
    //printf("return(x) = ");    fq_nmod_poly_print_pretty(result, "X", ctx); printf("\n\n");  // sqrt(Qx) ���
    
    fq_nmod_poly_t calre;                               // ��갪�� ���� ���� ������ ���ϱ� ����.  ��������� result�� ������ Q(x)�� �������� Ȯ��
    
    fq_nmod_poly_init(calre, ctx);
    fq_nmod_poly_mulmod(calre, result, result, g, ctx); //result�� ������ mod g(x)�� ���� ����Ͽ� calre�� ����.
    
    //printf("cal       = ");    fq_nmod_poly_print_pretty(calre, "X", ctx);  printf("\n");  //ri

    int a,b;
    
    a = fq_nmod_poly_equal(Qx,calre,ctx);               // Qx == calre �� ��� 1���, �ƴϸ� 0���

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
    
    //printf("time algori : %fs\n", t_precomp);     // �˰��� ���ð�

    //================== ������ ���� ================================
    // fqm[X]/g(X)�� ���� f(X)�� ���� f(X)^(2^(mt-1)) mod g(X)�� sqrt(f(X))�� ���� �̿�.
    fq_nmod_poly_t expone;          
    fq_nmod_poly_init(expone,ctx);


    c0 = clock();                       // time start
    fq_nmod_poly_powmod_fmpz_binexp(expone, Qx, twomt_1, g, ctx);     // Qx^(2^(mt-1)) mod g(X) 
    //printf("Sq(x) = ");     fq_nmod_poly_print_pretty(h,"X",ctx);    printf("\n");  //Sq Ȯ��
    c1 = clock();                       // time end
    
    
    t_modexp += (double) (c1 - c0) / CLOCKS_PER_SEC;

    b = fq_nmod_poly_equal(result,expone,ctx);          // ���� �˰���� ���� ����� �������� Ȯ��. 
    
    if(b==FALSE) //��� Ȯ�ο�
    {
        printf("\n FALSE \n");
        break;
    }
    else
        printf(" TRUE sqrt(Q(x)) == return \n");   

    //printf("time expone : %fs\n", t_modexp); // ������ ���ð�
    //printf("\n\n");
    
    //====================== �޸� ���� ============================
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
    printf("time algori : %fs\n", t_precomp/n_iter);     // �˰��� ���ð�
    printf("time expone : %fs\n", t_modexp/n_iter);     // ������ ���ð�


    return EXIT_SUCCESS;
}