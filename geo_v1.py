from math import *
import numpy as np

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a # splaszczenie
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # mimowsrod eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # mimosrod^2 eccentricity**2

    def deg2dms(self, dd):
        deg = int(np.trunc(dd))
        mnt = int(np.trunc((dd-deg) * 60))
        sec = ((dd-deg) * 60 - mnt) * 60
        dms = [deg, abs(mnt), abs(sec)]
        #print(str(deg)+chr(176)+"%0.2d" % abs(mnt)+'\''+"%08.5f" % abs(sec)+'\"')
        return dms
    
    def Np(self,lat):
        """
        """
        N = self.a/sqrt(1-self.ecc2*sin(lat))
        return N
    def Mp(self,lat):
        """
        """
        M = (self.a*(1-self.ecc2))/sqrt((1-self.ecc2*sin(lat)**2)**3)
        return M
    
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")
            
        def plh2xyz(self,P,L,H):
            """
            """
            N = self.Np(P)
            X = (N+H)*cos(P)*cos(L)
            Y = (N+H)*cos(P)*sin(L)
            Z = (N*(1 -self.ecc2 )+H)*sin(P)
            return X,Y,Z
        
        def Rneu(P,L):
            """
            """
            R = np.array([[-sin(P)*cos(L) , -sin(L), cos(P)*cos(L)]
                         ,[-sin(P)*sin(L) , cos(L) , cos(P)*sin(L)]
                         ,[cos(fi)        , 0      ,    sin(P)  ]])
            return R 


        def xyz2neu(P,L,dX): #dX is an array!
            """
            """
            R = Rneu(P,L)
            RT = np.transpose(R)
            dx = np.dot(RT,dX)
            return dx #array !
        
        
        def sigma(self,P):
            """
            """
            A0 = 1-(self.ecc2/4)-((3*(self.ecc2**2))/64)-((5*(self.ecc2**3))/256);
            A2 = (3/8)*(self.ecc2+(self.ecc2**2/4)+((15*(self.ecc2**3))/128));
            A4 = (15/256)*((self.ecc2**2)+((3*(self.ecc2**3))/4));
            A6 = (35*(self.ecc2**3))/3072;
            si= self.a*(A0*P-A2*sin(2*P)+A4*sin(4*P)-A6*sin(6*P));
            return si
        
        
        def pl2xy(self,P,L,L0):
            """
            """
            b2 = (self.a**2)*(1-self.ecc2);
            epr2 = ((self.a**2)-(b2))/(b2);
            t = tan(P);
            n2 = epr2*((cos(P))**2);
            
            N = self.Np(P)
            si = self.sigma(P)
            dl = L - L0
            
            xgk = si+ ((dl**2)/2)*N*sin(P)*cos(P)*(1+((dl**2)/12)*((cos(P))**2)*(5-t**2+9*n2+4*(n2**2))+((dl**4)/360)*((cos(P))**4)*(61-58*(t**2)+t**4+270*n2-330*n2*(t**2)))
            ygk = dl*N*cos(P)*(1+((dl**2)/6)*((cos(P))**2)*(1-t**2+n2)+((dl**4)/120)*((cos(P))**4)*(5-18*(t**2)+t**4+14*n2-58*n2*(t**2)))
            
            return xgk,ygk
        
        
        def f1(self,xgk):
            """
            """
            A0 = 1-(self.ecc2/4)-((3*(self.ecc2**2))/64)-((5*(self.ecc2**3))/256) 
            f = xgk/(self.a*A0) 
            
            while 1:
                si = self.sigma(f)
                fpomoc=f
                f = f + (xgk-si)/(self.a*A0)
                if abs(f-fpomoc)<(0.000001/206265):
                    break
            
            return f
        
        def xy2pl(self,xgk,ygk,L0):
            """
            """
            fp = self.f1(xgk)
            N = self.Np(fp)
            M = self.Mp(fp)

            t1 = tan(fp)
            b2 = (self.a**2)*(1-self.ecc2)
            epr2 = ((self.a**2)-(b2))/(b2)
            n12 = epr2*((cos(fp))**2)

            p= fp-((ygk**2)*t1/(2*M*N))*(1-((ygk**2)/(12*(N**2)))*(5+3*(t1**2)+n12-9*n12*(t1**2)-4*(n12**2))+((ygk**4)/(360*(N**4)))*(61+90*(t1**2)+45*(t1**4)));
            l= L0+(ygk/(N*cos(fp)))*(1-((ygk**2)/(6*(N**2)))*(1+2*(t1**2)+n12)+((ygk**4)/(120*(N**4)))*(5+28*(t1**2)+24*(t1**4)+6*n12+8*n12*(t1**2)));
            
            return p,l
        
        def u1992(xgk,ygk):
            """
            """
            x = xgk * 0.9993 - 5300000
            y = ygk * 0.9993 + 500000
            return x,y
        
        def u2000(xgk,ygk,L0): #L0 15st,18st,21st,24st
            """
            """
            x = xgk * 0.999923
            y = ygk * 0.999923 + (L0*180/pi/3) * 1000000 + 500000
            return x,y
        

if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    phi, lam, h = geo.xyz2plh(X, Y, Z)
    print(phi, lam, h)
    #phi, lam, h = geo.xyz2plh2(X, Y, Z)
    #print(phi, lam, h)
        
    
