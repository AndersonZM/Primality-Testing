//Anderson Zudio
//AC em 1/7/2014

#include<stdio.h>

int main()
{
    long int testcase, n, a, b, c,sum,i,j;
    while(scanf("%ld",&testcase)==1)
    {
        for(i=1;i<=testcase;i++)
        {
            scanf("%ld",&n);
            sum=0;
            for(j=1;j<=n;j++)
            {
                scanf("%ld%ld%ld",&a,&b,&c);
                sum+=(a*c);
            }
            printf("%ld\n",sum);
        }
    }
    return 0;

}
