import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.util.*;

public class VAAS_Sim implements PlugIn {
        double disFactor = 0.7;//distance factor between dots, 1*psf_width/2
        double AU = 0.1;//pinhole size
        public void run(String arg) {
            
            int width =1024;
            int height = 1024;
            double pixelSize = 10.0; // nm per pixelSize
            double lmda = 550; //nm
            double NA = 1.49;
            double psf_width = 1.22*lmda/NA;
            double I = 10; // initial intensity

            double[] image = new double[width*height];
            double[] vaas  = new double[width*height];

            int margin = (int)(psf_width/2.0/pixelSize+0.5);
            int dx = (int)(psf_width/2*disFactor/pixelSize + 0.5);
            int dy = (int)(psf_width/2*disFactor/pixelSize + 0.5);

            int x_num_pnts = 4;
            int y_num_pnts = 4;

            Point[] pnts = new Point[x_num_pnts*y_num_pnts];

            for (int j=0;j<y_num_pnts;j++)
                for (int i=0;i<x_num_pnts;i++) {
                    int idx = j*x_num_pnts+i;
                    pnts[idx] = new Point(100+i*dx,100+j*dy);
                }

            for (int n=0;n<x_num_pnts*y_num_pnts;n++) {
                    
                IJ.showProgress(n, x_num_pnts*y_num_pnts);
                
                int x0 = pnts[n].x;
                int y0 = pnts[n].y;
                for (int y = y0-margin;y<y0+margin;y++)
                    for (int x=x0-margin; x<x0+margin;x++) {
                        //printf("%d,%d\n",x,y);
                        int idx = y*width+x;
                        double I0 = I*Airy(Math.sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))/margin*3.83);
                        for (int j =y-margin;j<y+margin;j++)
                            for (int i=x-margin;i<x+margin;i++) {
                                double detRadius2 = (i-x)*(i-x)+(j-y)*(j-y);// R square. Avoid sqrt to save time.
                                double illumRadius2 = (i-x0)*(i-x0)+(j-y0)*(j-y0);
                                double Ie;//emission intensity
                                if(illumRadius2 < margin*margin)
                                    Ie = I0*Airy(Math.sqrt(illumRadius2)/margin*3.83);
                                else
                                    Ie = 0;
                                double pinhole = margin*AU;
                                if (detRadius2 > pinhole*pinhole)//use phinhole square to avoid sqrt
                                    vaas[idx]+=Ie;
                                else
                                    image[idx]+=Ie;
                            }
                    }
            }
            
            ImagePlus imp = IJ.createImage("Untitled", "16-bit black", width, height, 5);
            IJ.run(imp, "Set Scale...", "distance=1 known="+pixelSize+" unit=nm");
            imp.show();
            double maxVaas = findMax(vaas);
            double maxImage = findMax(image);
            for(int page=1;page<=5;page++){
                imp.setSlice(page);
                ImageProcessor ip = imp.getProcessor();
                for(int j=0;j<height;j++)
                    for(int i=0;i<width;i++){
                            int idx = i+j*width;
                            if(page==1)// image from pinhole
                                ip.putPixel(i,j,(int)(image[idx]));
                            if(page==2)// image from vaas
                                ip.putPixel(i,j,(int)(vaas[idx]));
                            if(page==3){//image devided by vaas image
                                double f = 1.0-vaas[idx]/maxVaas;
                                ip.putPixel(i,j,(int)(f*image[idx]));
                            }
                            if(page==4){//the supper resolution, 0.1AU*1AU
                                double f = image[idx]/maxImage;
                                ip.putPixel(i,j,(int)(f*(image[idx]+vaas[idx])));
                            }
                            if(page==5){//integration from vaas ring, does not work well
                                int vaasRadius = (int)(margin/2.0);
                                int IRing = 0;
                                for(int degree=0;degree<360;degree+=18){
                                        double a = degree/180*3.14;
                                        int ringX = i+(int)(vaasRadius*Math.cos(a));
                                        int ringY = j+(int)(vaasRadius*Math.sin(a));
                                        int ringIdx = ringX+ringY*width;
                                        if(ringIdx > 0 && ringIdx < width*height)
                                                IRing += vaas[ringX+ringY*width]; 
                                }
                                ip.putPixel(i,j,IRing);
                            }
                    }
            }

            
        }
        /*********************************************/
        double findMax(double[] array){
                double max = -1.0/0 ;
                for(int i=0; i<array.length;i++){
                        if (max < array[i]) max = array[i];
                }
                return max;
        }
        /*********************************************/
        double Airy( double x ) {
            if(x==0) return 1.0;
            double f0= 2*J1(x)/x;
            return f0*f0;
        }
        /*********************************************/
        public static double J1(final double x)
        {
            double ax;

            if ((ax = Math.abs(x)) < 8.0)
            {
                final double y = x * x;
                final double ans1 = x *
                                    (72362614232.0 + y *
                                     (-7895059235.0 + y *
                                      (242396853.1 + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
                final double ans2 = 144725228442.0 + y *
                                    (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
                return ans1 / ans2;
            }

            final double z = 8.0 / ax;
            final double xx = ax - 2.356194491;
            final double y = z * z;

            final double ans1 = 1.0 + y *
                                (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
            final double ans2 = 0.04687499995 + y *
                                (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
            final double ans = Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * ans1 - z * Math.sin(xx) * ans2);
            return (x < 0.0) ? -ans : ans;
        }
}





