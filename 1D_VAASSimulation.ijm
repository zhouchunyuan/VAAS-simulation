width = 512;
height = 512;

pixelSize = 3; // nm per pixel
lmda = 550; // nm
NA = 1.0;
PSF_WIDTH = 0.61*lmda/NA / pixelSize;
sigma = PSF_WIDTH/6;

I = 65535;// initail intensity
title = "confocal simulator";
/***********************************
 if we set pinhole = 1.2 AU, 
 the two dots will appear to
 merge at 0.4 AU distance

 if we set pinhole = 0.1 AU,
 two peaks can be seen even we
 we move them as close as 0.28 AU
***********************************/
pinholeSizeFactor = 0.5;//0.1 - 1
distanceFactor =0.26;// 0.25 - 0.35

distance = floor(PSF_WIDTH*distanceFactor);//distance between two dots
offset = 0;// center offset of these dots

// if pinhole larger,attenuation should be smaller
// and if distance closer, also attenuation should be smaller 
attenuationFactor = 0.008/pinholeSizeFactor*distanceFactor;

newImage(title, "16-bit black", 512, 512, 1);
//cmdstr = "distance=1 known="+pixelSize+"unit=nm";
run("Set Scale...", "distance=1 known="+pixelSize+" unit=nm");

illum = newArray(width);
fluo = newArray(width);
detect = newArray(width);
detect_shift= newArray(width);
vaas_out = newArray(width);

p1 = floor(width/2)-floor(distance/2) + offset;
p2 = floor(width/2)+floor(distance/2) + offset;

for(illum_center = 0; illum_center< width; illum_center++){ 
    // create illumination array
    for(i=0;i<width;i++){
        x = (i-illum_center)/PSF_WIDTH*2*1.22*PI;
        //illum = I*exp(-(i-illum_center)*(i-illum_center)/(2*sigma*sigma));
        illum[i] = I*Airy(x);

        setPixel(i,height/2,illum[i]);
    }
    // create fluorescent array
    for(i=0;i<width;i++){
        x1 = (i-p1)/PSF_WIDTH*2*1.22*PI;
        x2 = (i-p2)/PSF_WIDTH*2*1.22*PI;
        I1 = illum[floor(p1)];// illumination intensity for p1
        I2 = illum[floor(p2)];// illumination intensity for p2
        fluo[i] = I1*Airy(x1)+I2*Airy(x2);
        setPixel(i,height/2+10,fluo[i]);
    }
    // create center detector array
    AU = PSF_WIDTH*pinholeSizeFactor ;
    detect_pos = illum_center;
    signal = 0;
    for(i=floor(detect_pos-AU/2);i<floor(detect_pos+AU/2);i++){
        if(i<0)
            signal +=fluo[0];
        else if(i>=width)
            signal +=fluo[width-1];
        else
            signal +=fluo[i];
    }
    detect[detect_pos] = signal*attenuationFactor ;
    setPixel(detect_pos,height/2+50,detect[detect_pos]);

    // create shifted detector array
    /*********************************
      here we simulate a PMT located
      an bias-distance away from the center.
      We can see, a shifted PMT can increase 
      resolution, but decrease the intensity.
     *********************************/
    bias_distance = PSF_WIDTH*0.5;
    detect_pos = illum_center + bias_distance;
    signal = 0;
    for(i=floor(detect_pos-AU/2);i<floor(detect_pos+AU/2);i++){
        if(i<0)
            signal +=fluo[0];
        else if(i>=width)
            signal +=fluo[width-1];
        else
            signal +=fluo[i];
    }
    detect_shift[illum_center] = signal ;
    setPixel(illum_center,height/2+100,detect_shift[illum_center]*attenuationFactor); 

    // create VAAS condition
    detect_pos = illum_center;
    signal_out = 0;
    for(i=floor(detect_pos-PSF_WIDTH/2);i<floor(detect_pos+PSF_WIDTH/2);i++){
        if((i< floor(detect_pos-AU/2)) || (i>floor(detect_pos+AU/2))){
            if(i<0)
                signal_out +=fluo[0];
            else if(i>=width)
                signal_out +=fluo[width-1];
            else
                signal_out +=fluo[i];
        }
    }
    vaas_out[detect_pos] = signal_out ;
    setPixel(illum_center,height/2+150,vaas_out[illum_center]*attenuationFactor );
   
    updateDisplay();
    //wait(10);
}
/************
profileAt = height/2+100;
makeLine(0,profileAt ,width,profileAt );
run("Plot Profile");
selectWindow(title);
profileAt = height/2+150;
makeLine(0,profileAt ,width,profileAt );
run("Plot Profile");
*************/

y_axis = newArray(width);
x_axis = newArray(width);
Array.getStatistics(vaas_out, min_vaas, max_vaas, mean, stdDev);
Array.getStatistics(detect, min_detect, max_detect, mean, stdDev);
for(i=0;i<width;i++){
    x_axis[i]=i;
    if(detect[i] > max_detect/100)
        y_axis[i] = detect[i]/(vaas_out[i]/max_vaas);
    else
        y_axis[i] = detect[i];
}

Plot.create("Ratio", "Distance","Ratio", x_axis,y_axis);
Plot.show();

for(i=0;i<width;i++){
    x_axis[i]=i;
    y_axis[i] = detect_shift[i];
}

Plot.create("shifted detector", "Distance","I", x_axis,y_axis);
Plot.show();
/******************************************************
# Since there is no Bessel function in Imagej macro,
# I have to use j1 = sin(x)/(x*x)-cos(x)/x to simulate
# Bessel first kind first order J1
# Scale factor of 1.14 in x, and 2.25 in y is found by
# python 3.6:
# Comparison vs gaussian is also done. With sigma of 3.83/3
# The gaussian is very similar to Airy disk

import numpy as np
import pylab as py
import scipy.special as sp

def airy(x):
    f = 1.14
    a = np.sin(x*f)/((x*f)**2)-np.cos(x*f)/(x*f)
    return 2.25*((2*a/(x*f))**2)

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

x = np.linspace(-10, 10, 2000)
py.plot(x, (2*sp.j1(x)/x)**2,'r--',
        x, airy(x),'b--',
        x, gaussian(x,0,3.83/3),'g--')

py.xlim((-10, 10))
py.ylim((-0.5, 1.1))
py.legend(('$(2\mathcal{J}_1(x)/x)^2$',
           r'$\mathrm{2.25}(\frac{sin(x)/x^2-cos(x)/x}{x})^2$',
           '$gaussian$'))
py.xlabel('$x$')
py.ylabel('Intensity')
py.grid(True)
                                     
py.show()
*******************************************************/
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

