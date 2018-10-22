N = 1000;
y1 = newArray(N);
y2 = newArray(N);
sum = newArray(N);
x = newArray(N);
range = 30;
step = range/N;
d = 3.83/2;
for(i=0;i<N;i++){
r = i*step-range/2;
y1[i] = Airy(r-d);
y2[i] = Airy(r+d);
sum[i] = y1[i]+y2[i];
x[i] = r;
}

Plot.create("test","x","y",x,sum);
Plot.add("line",x,y2);
Plot.add("line",x,y1);
Plot.show();

/****************************************************/
function AiryByExp(x){
    s = 3.83/3;
    return exp(-x*x/(2*s*s));
}
function Airy(x){
    if(x==0)return 1;
    j = J1(x);
    return 4*j*j/x/x;
}
function J1(x)
{
            ax = abs(x);

            if (ax < 8.0)
            {
                y = x * x;
                ans1 = x *
                                    (72362614232.0 + y *
                                     (-7895059235.0 + y *
                                      (242396853.1 + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
                ans2 = 144725228442.0 + y *
                                    (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
                return ans1 / ans2;
            }

            z = 8.0 / ax;
            xx = ax - 2.356194491;
            y = z * z;

            ans1 = 1.0 + y *
                                (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
            ans2 = 0.04687499995 + y *
                                (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
            ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
            if(x<0.0) return -ans;
            else return ans;
}
