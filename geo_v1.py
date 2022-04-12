from math import *
import numpy as np

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            L0_2000 - południk osiowy dla układu PL-2000 [rad]
            L0_1992 - południk osiowy dla układu PL-1992 [rad]
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        self.L0_2000 = 21*pi/180
        self.L0_1992 = 19*pi/180
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
        """
        Parameters
        ----------
        dd : float - stopnie w liczbie dziesiętnej [stopnie]

        Returns
        -------
        dms : tuple - stopnie, minuty, sekundy

        """
        deg = int(np.trunc(dd))
        mnt = int(np.trunc((dd-deg) * 60))
        sec = ((dd-deg) * 60 - mnt) * 60
        dms = deg, abs(mnt), abs(sec)
        #print(str(deg)+chr(176)+"%0.2d" % abs(mnt)+'\''+"%08.5f" % abs(sec)+'\"')
        return dms
    
    def Np(self,lat):
        """
        Max. promień w I wertykale w kierunku głownym.
        Parametry
        ----------
        lat : float - szerokosc geodezyjna [radians]
        
        Returns
        -------
        N : float - promien w I wertykale [m]
        """
        N = self.a/sqrt(1-self.ecc2*sin(lat)**2)
        return N
    def Mp(self,lat):
        """
        Min. promień w kwadracie I mimosrodu w kierunku głównym.
        Parametry
        ---------
        lat : float - szerokosc geodezyjna [radians]
        
        Returns
        -------
        M : float - promien w kwadracie I mimosrodu
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
            
    def plh2xyz(self,P,L,H, inp = 'dec_degree'):
        """
        Transformacja odwrotna do algorytmu Hirvonena- transformacja współrzędnych geodezyjnych:
        szerokosc, dlugosc i wysokosc elipsoidalna na współrzędne ortokartezjańskie X,Y,Z
        Parametres:
        ----------
        P,L,H - float - szerokosc, dlugosc oraz wysokosc elipsoidalna [stopnie dziesietne i metry]
        inp - [str] , optional: dec_degree, radians
        Returns:
        --------
        X,Y,Z - float - wspolrzedne ortokartezjańskie [m]
        
        """
        if inp == 'dec_degree':
            P = P*pi/180
            L = L*pi/180
            N = self.Np(P)
            X = (N+H)*cos(P)*cos(L)
            Y = (N+H)*cos(P)*sin(L)
            Z = (N*(1 -self.ecc2 )+H)*sin(P)
    
        elif inp == 'radians':
            N = self.Np(P)
            X = (N+H)*cos(P)*cos(L)
            Y = (N+H)*cos(P)*sin(L)
            Z = (N*(1 -self.ecc2 )+H)*sin(P)
        return X,Y,Z
        
    def Rneu(self,P,L):
        """
        Macierz obrotu R
        Parametres:
        ----------
        P,L - float - szerokosc i dlugosc geodezyjna [radians]
        
        Returns:
        --------
        R - array - macierz obrotu
        
        """
        R = np.array([[-sin(P)*cos(L) , -sin(L), cos(P)*cos(L)]
                      ,[-sin(P)*sin(L) , cos(L) , cos(P)*sin(L)]
                      ,[cos(P)        , 0      ,    sin(P)  ]])
        return R 


    def xyz2neu(self,P,L,dX,inp = 'dec_degree'): #dX is an array!
        """
        Macierz delta_N,delta_E,delta_U. Posłużenie się funkcją Rneu.
        Parametres:
        ----------
        P,L - float - współrzędne elipsoidalne punktu referencyjnego dla układu topocentrycznego [stopnie dziesietne]
        dX - array - różnice między współrzędnymi kolejnych punktów, a współrzędnymi punktu referencyjnego w układzie ortokartezjańskim [m]
        
        Returns:
        --------
        dx - array - macierz delta_N,delta_E,delta_U [m]
        """
        if inp == 'dec_degree':
            P = P*pi/180
            L = L*pi/180
        R = self.Rneu(P,L)
        RT = np.transpose(R)
        dx = RT @ dX
        return dx #array !
        
        
    def sigma(self,P):
        """
        Wartosc sigmy
        Paramtetres:
        -----------
        P - float - szerokosc elipsoidalna [radians]
        
        Returns:
        --------
        si - float - sigma
        """
        A0 = 1-(self.ecc2/4)-((3*(self.ecc2**2))/64)-((5*(self.ecc2**3))/256);
        A2 = (3/8)*(self.ecc2+(self.ecc2**2/4)+((15*(self.ecc2**3))/128));
        A4 = (15/256)*((self.ecc2**2)+((3*(self.ecc2**3))/4));
        A6 = (35*(self.ecc2**3))/3072;
        si= self.a*(A0*P-A2*sin(2*P)+A4*sin(4*P)-A6*sin(6*P));
        return si
        
        
    def pl2xy(self,P,L,L0, inp = 'dec_degree'):
        """
        Transformacja wspolrzednych geodezyjnych na wspolrzedne w układzie Gaussa-Krugera
        Parametres:
        P,L - float - szerokosc oraz dlugosc geodezyjna [stopnie dziesietne]
        L0 - float - poludnik osiowy [radians]
        
        Returns:
        xgk,ygk - float - wspolrzedne w ukladzie Gaussa-Krugera [m]
        """
        if inp == 'dec_degree':
            P=P*pi/180
            L=L*pi/180
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
        Funkcja pomocnicza do obliczenia wspolrzednych geodezyjnych z dokladnoscia do
        0.000001 sekundy
        Parametres:
        -----------
        xgk - wspolrzedna X w ukladzie Gaussa-Krugera [m]
        
        Returns:
        --------
        f - wartosc pomocnicza do szerokosci geodezyjnej
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
        
    def xy2pl(self,xgk,ygk,L0, output = 'dec_degree'):
        """
        Transformacja wspolrzednych w ukladzie Gaussa-Krugera na wspolrzedne geodezyjne. Wykorzystana
        funkcja f1
        Parametres:
        -----------
        xgk,ygk - float - wspolrzedne w ukladzie Gaussa-Krugera [m]
        L0 - float - poludnik osiowy [radians]
        
        Returns:
        --------
        p,l - float - wspolrzedne geodezyjne [stopnie dziesietne]
        output - optional, default
        dec_degree - stopnie dziesietne
        dms - degrees, minutes, seconds
        """
        fp = self.f1(xgk)
        N = self.Np(fp)
        M = self.Mp(fp)

        t1 = tan(fp)
        b2 = (self.a**2)*(1-self.ecc2)
        epr2 = ((self.a**2)-(b2))/(b2)
        n12 = epr2*((cos(fp))**2)

        p= fp-((ygk**2)*t1/(2*M*N))*(1-((ygk**2)/(12*(N**2)))*(5+3*(t1**2)+n12-9*n12*(t1**2)-4*(n12**2))+((ygk**4)/(360*(N**4)))*(61+90*(t1**2)+45*(t1**4)))
        l= L0+(ygk/(N*cos(fp)))*(1-((ygk**2)/(6*(N**2)))*(1+2*(t1**2)+n12)+((ygk**4)/(120*(N**4)))*(5+28*(t1**2)+24*(t1**4)+6*n12+8*n12*(t1**2)))
            
        if output == 'dec_degree':
            return degrees(p),degrees(l)
        elif output == 'dms':
            lat = self.deg2dms(degrees(p))
            lon = self.deg2dms(degrees(l))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}" 
                
                
        
    def u1992(self,xgk,ygk):
        """
        Transformacja współrzędnych w układzie Gaussa-Krugera do układu PL-1992
        Parametres:
        -----------
        xgk,ygk - float - współrzędne w układzie Gaussa-Krugera [m]
        
        Returns:
        --------
        x,y - float - współrzędne w układzie PL-1992 [m]
        """
        x = xgk * 0.9993 - 5300000
        y = ygk * 0.9993 + 500000
        return x,y
        
    def u2000(self,xgk,ygk,L0 = 21*pi/180): #L0 15st,18st,21st,24st
        """
        Transformacja współrzędnych w układzie Gaussa-Krugera do układu PL-2000
        Parametres:
        -----------
        xgk,ygk - float - współrzędne w układzie Gaussa-Krugera [m]
        L0 - float - poludnik osiowy [radians]
        
        Returns:
        --------
        x,y - float - współrzędne w układzie PL-2000 [m]
        """
        x = xgk * 0.999923
        y = ygk * 0.999923 + (L0*180/pi/3) * 1000000 + 500000
        return x,y
    
    def s_azimuth_elevation(self,dN,dE,dU):
        """
        Obliczenie odległosci 3D, azymutu oraz kąta elewacji dla punktu kolejnego oraz punktu
        odniesienia w ukladzie topocentrycznym
        Parametres:
        -----------
        dN,dE,dU - float - współrzędne w układzie topocentrycznym [m]
        
        Returns:
        --------
        s - float - odleglosc 3D [m]
        Az,z - float - azymut oraz kąt elewacji [stopnie dziesietne]
        """
        s = sqrt(dN**2 + dE**2 + dU**2)
        Az = atan2(dE, dN)
        z = acos(dU/s)
        
        if Az > 0 and z > 0:
            Az = (Az)*180/pi
            z = (z)*180/pi
        elif Az < 0 and z < 0:
            Az = (Az+2*pi)*180/pi
            z = (z+2*pi)*180/pi
        elif Az < 0 and z > 0:
            Az = (Az+2*pi)*180/pi
            z = z*180/pi
        elif Az > 0 and z < 0:
            Az = Az*180/pi
            z = (z+2*pi)*180/pi            
        return s,Az,z
