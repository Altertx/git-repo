from math import *
class Transformacje:
    def __init__(self, model: str = "wgs84"):
    #https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
             raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flattening - self.flattening ** 2)
    def xyz2plh(self, X, Y, Z):
        #wkleic sobie hirvonena

        return phi,lam,h
        
if __name__ == "__main__":
    #utworzenie obiektu
    geo = Transformacje(model = "grs80")
    #dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170