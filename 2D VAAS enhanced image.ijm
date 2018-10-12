/******** basic parameters *********/
width = 1024;
height = 1024;

pixelSize = 10; // nm per pixel
lmda = 550; // nm
NA = 1.0;
PSF_WIDTH = 1.22*lmda/NA / pixelSize;

I = 100; // initial intensity
/************************************/

image = newArray(width*height);
vaas = newArray(width*height);

margin = floor(PSF_WIDTH/2);

dx = PSF_WIDTH*0.35;
dy = PSF_WIDTH*0.35;

x_num_pnts = 4;
y_num_pnts = 4;
pntsX = newArray(x_num_pnts*y_num_pnts);
pntsY = newArray(x_num_pnts*y_num_pnts);
for(j=0;j<y_num_pnts;j++){
    for(i=0;i<x_num_pnts;i++){
        idx = j*x_num_pnts+i;
        pntsX[idx] = 100+i*dx;
        pntsY[idx] = 100+j*dx;

    }
}

pixArraySize = margin*margin;

for(n=0;n<lengthOf(pntsX);n++){
    x0 = margin + floor(pntsX[n]);
    y0 = margin + floor(pntsY[n]);
    
    showProgress(n , lengthOf(pntsX));
    
    for(y =y0-margin; y <y0+margin; y++){
        for(x =x0-margin; x <x0+margin; x++){
            pxIndex = y*width+x;
            /********* this is the ideal confocal **************
            //* I0 is the illumination intensity on this dot 
            I0 = I*exp(-(pow(x-x0,2)+pow(y-y0,2))/(2*pow(margin/3,2)));
            //* I0 = I*Airy(sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))/margin*3.83);
            Ie = I0*I0/I;// Ie is the fluorescent intensity on (x,y)
            image[y*width+x] += Ie;
            ***************************************************/
            I0 = I*Airy(sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))/margin*3.83);

            for(j =y-margin; j <y+margin; j++){

                for(i =x-margin; i <x+margin; i++){
                    R = sqrt((i-x)*(i-x)+(j-y)*(j-y));
                    Ie = I0*Airy(sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0))/margin*3.83);
                    if(R>margin*0.2)
                        vaas[pxIndex ] += Ie;
                    else
                        image[pxIndex ] += Ie;

                }

            }

        }
    }
}

newImage("test", "16-bit black", width, height, 4);

setSlice(1);
for(j=0;j<height;j++){
    for(i=0;i<width;i++){
        setPixel(i,j,image[i+j*width]);
    }
}
setSlice(2);
for(j=0;j<height;j++){
    for(i=0;i<width;i++){
        setPixel(i,j,vaas[i+j*width]);
    }
}

Array.getStatistics(vaas, min, max, mean, stdDev);

setSlice(3);//division methord
for(j=0;j<height;j++){
    for(i=0;i<width;i++){
        idx = i+j*width;
        f = 1 - vaas[idx]/max;
        setPixel(i,j,f*image[idx]);
    }
}

Array.getStatistics(image, min, max, mean, stdDev);

setSlice(4);//multiplication methord
for(j=0;j<height;j++){
    for(i=0;i<width;i++){
        idx = i+j*width;
        f = image[idx]/max;
        setPixel(i,j,f*(image[idx]+vaas[idx]));
    }
}
resetMinAndMax();
/***********************/
function Airy(x){
    factor = 1.14;
    m = 2.25;
    fx = factor*x;
    if(fx==0){
        ret = 1;
    }else{
        j1 = sin(fx)/(fx*fx)-cos(fx)/fx;
        ret = m*(2*j1/fx)*(2*j1/fx);
    }
    return ret;
}
