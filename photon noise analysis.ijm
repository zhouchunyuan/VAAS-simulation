N= 1024;
photons = 1;
a = newArray(N*N);
b = newArray(N*N);
c = newArray(N*N);

for(j=0;j<N;j++){
    for(i=0;i<N;i++){
        idx = j*N+i;
        signal = photons;
        noise = sqrt(signal)*(1-random());
        a[idx] = signal+noise;

        signal = 0.8*photons;
        noise = sqrt(signal)*(1-random());
        b[idx] = signal+noise;

        c[idx] = a[idx]*b[idx];
    }
}

Array.getStatistics(a, min, max, mean, stdDev);
print("array a mean,std,error :(",mean,"),(",stdDev,")",stdDev/mean);

Array.getStatistics(b, min, max, mean, stdDev);
print("array b mean,std,error :(",mean,"),(",stdDev,")",stdDev/mean);

Array.getStatistics(c, min, max, mean, stdDev);
print("array c mean,std,error :(",mean,"),(",stdDev,")",stdDev/mean);
