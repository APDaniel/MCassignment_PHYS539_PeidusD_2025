using MCassignment_PHYS539_PeidusD_2025.Containers;
using System.Numerics;

namespace MCassignment_PHYS539_PeidusD_2025.Helpers
{
    public static class PhysicsData
    {
        public class CrossSections
        {
            public double Photoelectric;
            public double Compton;
            public double Pair;
            public double Total => Photoelectric + Compton + Pair;
        }
        
        public static CrossSections GetCrossSections(double energy)
        {
            if (Math.Abs(energy - 2.0) < 0.1)
               return new CrossSections { Photoelectric = 1.063e-6, Compton = 4.901e-2, Pair = 3.908e-4 };
            if (Math.Abs(energy - 6.0) < 0.1)
                return new CrossSections { Photoelectric = 2.483e-7, Compton = 2.455e-2, Pair = 3.155e-3 };
            if (Math.Abs(energy - 10.0) < 0.1)
                return new CrossSections { Photoelectric = 1.386e-7, Compton = 1.710e-2, Pair = 5.090e-3 };

            throw new ArgumentException("Unsupported photon energy");
        }
        public static double GetTotalLinearAttenuationCoefficient(double energy)
        {
            if (Math.Abs(energy - 2.0) < 0.1)
                return 0.049;
            if (Math.Abs(energy - 6.0) < 0.1)
                return 0.029;
            if (Math.Abs(energy - 10.0) < 0.1)
                return 0.022;

            throw new ArgumentException("Unsupported attenuation coefficient");
        }
        public static string SampleInteractionType(double energy)
        {
            var cs = GetCrossSections(energy);
            double totalCs = cs.Total;
            double peProb = cs.Photoelectric/ cs.Total;
            double pairProb = cs.Pair/ cs.Total;
            double comptonProb = cs.Compton/ cs.Total;
            var random = new Random();
            var R = random.NextDouble();
            if (R <= peProb) return "photoelectric";
            if (R > peProb && R <= pairProb) return "pair";
            if (R > pairProb && R <= 1) return "compton";
            else throw new ArgumentException("Check logic for interaction probabilities");
        }
        public static (double electronEnergy, Vector3 electronDirection, double positronEnergy, Vector3 positronDirection) GetElectronFromInteraction(
            string interaction, double photonEnergy, Vector3 photonDir, 
            StatisticsToExport statistics, double maxKleinNichProb)
        {
            double theta = 0;
            double phi = 0;
            double theta_p = 0;
            double phi_p = 0;
            Vector3 eDirection = photonDir;
            Vector3 pDirection = new Vector3(0,0,0);
            double eEnergy = 0;
            double pEnergy = 0;
            switch (interaction)
            {
                case "photoelectric":
                    eEnergy = Math.Round(photonEnergy,5);
                    eDirection = photonDir;
                    break;
                case "compton":
                    // Step 1: Sample realistic scattered photon direction
                    Vector3 scatteredPhotonDir = PhysicsData.SampleScatteredComptonPhotonDirection(photonEnergy, photonDir, maxKleinNichProb,statistics);

                    // Step 2: Calculate angle between original and scattered direction
                    double anglePhoton = Vector3AngleBetween(photonDir, scatteredPhotonDir);

                    double scatteredEnergy = photonEnergy / (1 + (photonEnergy / 0.511) * (1 - Math.Cos(anglePhoton)));

                    // Step 3: Compute electron energy
                    eEnergy = photonEnergy - scatteredEnergy;
                    statistics.SaveComptonEnergy(eEnergy);
                    // Step 4: Estimate electron angle based on momentum conservation
                    double cotThetaElectron = (1+ photonEnergy / 0.511)*Math.Tan(anglePhoton*0.5);

                    double thetaElectron = Math.Atan(1/cotThetaElectron);

                    // Step 5: Sample random azimuthal angle for spread
                    double phiElectron = PhysicsData.SampleAzimutalAngle();
                    Vector3 localEDir = PhysicsData.SphericalToCartesian(thetaElectron, phiElectron);
                    eDirection = PhysicsData.RotateToDirection(localEDir, photonDir);
                    statistics.RecordComptonScatteringAngle(anglePhoton);
                    break;

                case "pair":
                    double T = photonEnergy - 1.022;
                    pEnergy = SamplePairProductionEnergy(photonEnergy);
                    eEnergy = T - pEnergy;
                    //theta = SamplePairProductionAngleTranslation(eEnergy);
                    theta = SamplePairProductionAngleRejection(eEnergy); 


                    phi = SampleAzimutalAngle();

                    //theta_p = SamplePairProductionAngleTranslation(pEnergy);
                    theta_p = SamplePairProductionAngleRejection(pEnergy);
                    phi_p = SampleAzimutalAngle();

                    eDirection = SphericalToCartesian(theta, phi);
                    pDirection = SphericalToCartesian(theta_p, phi_p);

                    statistics.RecordPairScatteringAngle_e(theta);
                    statistics.RecordPairScatteringAngle_p(theta_p);
                    break;
            }
            return (eEnergy, eDirection, pEnergy, pDirection);
        }
        public static double SampleAzimutalAngle()
        {
            var random = new Random();
            var R = random.NextDouble();
            return 2*Math.PI*R;
        }
        public static double ComptonScatteredPhotonEnergy(double energy,double phi)
        {
            double mec2 = 0.511;
            return energy / (1 + (energy / mec2) * (1 - Math.Cos(phi)));
        }
        public static double SamplePairProductionEnergy(double photonEnergy)
        {
            double output;
            double Tmax = (photonEnergy - 1.022)/2;
            while (true)
            {
                var random = new Random();
                var R = random.NextDouble();
                var random1 = new Random();
                var R1 = random1.NextDouble();

                double x = (R*Tmax)/(photonEnergy-2*0.511);
                double fx = Math.Log(1000 * x + 1);
                double y = R1 * Math.Log(1000* Tmax + 1);
                output = R * Tmax;
                if (y < fx) 
                    return output;
            }
        }
        public static Vector3 SampleScatteredComptonPhotonDirection(
            double photonEnergyMeV, Vector3 incidentDirection, 
            double maxKleinNichProb, StatisticsToExport statistics)
        {
            Random rng = new Random();
            double theta, phi;
            var electronRadius = Constants.electronRadius_m;

            // Rejection sampling using Klein-Nishina formula
            while (true)
            {
                theta = rng.NextDouble() * Math.PI; // phi_p in [0, pi]
                var hv = photonEnergyMeV;
                var hv_prime = hv / (1 + (hv / 0.511) * (1 - Math.Cos(theta)));
                var a = Math.PI * Math.Sin(theta) * Math.Pow(electronRadius, 2);
                var b = Math.Pow((hv_prime / hv), 2);
                var c = hv / hv_prime + hv_prime / hv - Math.Pow(Math.Sin(theta), 2);
                var cs = a * b * c;

                Random random = new Random();
                var R2 = random.NextDouble();
                double sampledMaxProb = R2*maxKleinNichProb; // Safe upper bound

                if (sampledMaxProb <= cs)
                    break;
            }

            // Azimuthal angle sampled uniformly
            phi = 2 * Math.PI * rng.NextDouble();

            // Local direction (relative to z-axis)
            Vector3 localDirection = SphericalToCartesian(theta, phi);
            
            // Rotate local direction to align with incident photon direction
            return RotateToDirection(localDirection, incidentDirection);
        }
        public static Vector3 RotateToDirection(Vector3 local, Vector3 targetDir)
        {
            Vector3 zAxis = new Vector3(0, 0, 1);
            if (Vector3.Distance(Vector3.Normalize(targetDir), zAxis) < 1e-6)
                return local;

            Vector3 axis = Vector3.Cross(zAxis, targetDir);
            float angle = (float)Math.Acos(Vector3.Dot(zAxis, Vector3.Normalize(targetDir)));
            Quaternion rotation = Quaternion.CreateFromAxisAngle(Vector3.Normalize(axis), angle);
            return Vector3.Transform(local, rotation);
        }
        public static double SamplePairProductionAngleRejection(double particleEnergy)
        {
            //Sample scattering angle randomly from 0 to PI
            var random = new Random();
            var R = random.NextDouble();
            R = R * Math.PI;

            double beta = Math.Sqrt(1 - 1 / Math.Pow((particleEnergy / 0.511) + 1, 2));
            double maxPairProdCs = PairProductionHelper.SamplePairProductionCrSct(particleEnergy);
            double alarama;
            double scatteringAngle=0.0;
            if (beta > 1 || beta < 0)
                alarama = 0.0;
            
            while (true)
            {
                var random1 = new Random();
                var R1 = random1.NextDouble();
                var phi = R1 * Math.PI;

                var thisPairProdCs = (Math.PI * Math.Sin(phi)) / Math.Pow((1 - beta * Math.Cos(phi)), 2);

                var random2 = new Random();
                var R2 = random2.NextDouble();
                R2 = R2 * maxPairProdCs;

                if (thisPairProdCs > R2)
                {
                    scatteringAngle = phi;
                    return scatteringAngle;
                }
            }
        }
        public static double SamplePairProductionAngleTranslation(double particleEnergy)
        {
            //Sample scattering angle randomly from 0 to PI
            var random = new Random();
            var R = random.NextDouble();
            R = R * Math.PI;

            double beta = Math.Sqrt(1 - 1 / Math.Pow((particleEnergy / 0.511) + 1, 2));

            var up = (1 + beta * R - R);
            var down = (1 - beta * R + (Math.Pow(beta, 2)) * R);
            double cosAngle = up / down;
            double angle = Math.Acos(cosAngle);

            return angle;
        }
        public static Vector3 SphericalToCartesian(double theta,double phi)
        {
            double sinTheta = Math.Sin(theta);
            return new Vector3(
                (float)(sinTheta*Math.Cos(phi)),
                (float)(sinTheta * Math.Sin(phi)),
                (float)(Math.Cos(theta)));
        }
        public static double SampleMultipleScatteringAngleRejection (double eMeV,double pathLength)
        {
            double thetaSquared = GetMeanSquaredScatteringAngleSquaredForWater(eMeV,pathLength);
            double maxProbability =maxCollisionScatteringProbability.GetMaxColScatteringProb(eMeV,pathLength);
            if (eMeV > 0)
            {
                while (true)
                {
                    var random = new Random(); 
                    var R1 = random.NextDouble();
                    var random1 = new Random(); 
                    var R2 = random1.NextDouble();

                    var a = 2 * R1 / Math.Pow(Math.PI, 2); 
                    a = a / thetaSquared;

                    var b = Math.Pow(R1, 2); 
                    b = b / thetaSquared; 
                    b = Math.Exp(-b);

                    var thisScatteringProbability = a * b;

                    var sampledMaxProb = R2 * maxProbability;

                    if (sampledMaxProb < thisScatteringProbability)
                        return R1;
                }
            }
            return 0.0;
        }
        public static double SampleMultipleScatteringAngleTranslation(double eMeV, double pathLength)
        {
            double thetaSquared = GetMeanSquaredScatteringAngleSquaredForWater(eMeV, pathLength);
            //double maxProbability = maxCollisionScatteringProbability.GetMaxColScatteringProb(eMeV, pathLength);
            if (eMeV > 0)
            {
                var random = new Random();
                var R1 = random.NextDouble();

                var log = Math.Log(1 - R1);
                var ratio = log * thetaSquared * (-1);
                var angle = Math.Sqrt(ratio);
                var random1 = new Random();
                var R2 = random1.NextDouble();
                if (R2 >= 0.5)
                    angle = angle * (-1);
                return angle;
            }
            return 0.0;
        }
        public static double GetMeanSquaredScatteringAngleSquaredForWater(double eMeV, double pathLength)
        {
            //Const for water
            double Z = Constants.Z;
            double A = Constants.A;
            double z = Constants.z;
            double m = Constants.m;
            double F = 0.9;
            var x = pathLength;

            var totalEnergy = eMeV + 0.511;
            var p = Math.Sqrt(Math.Pow(totalEnergy, 2) - Math.Pow(m, 2));
            var beta = p / totalEnergy;
            var Xc2 = 0.157 * z * ((Z * (Z + 1)) / A) * x / (Math.Pow(p, 2) * Math.Pow(beta, 2));
            double alpha = 0.00729927007299270072992700729927;

            var aa = (2.007e-5);
            double power = 2; power = power / 3;
            var bb = Math.Pow(Z, power);
            var cc = Math.Pow((Z * z * alpha / beta), 2);
            var Xa2 = aa * bb * (1 + 3.14 * cc) / Math.Pow(p, 2);

            var omega =Xc2/Xa2;
            var nu = 0.5 * omega / (1 - F);
            var a = 2*Xc2/(1+Math.Pow(F,2));
            var b = (1 + nu) / nu;
            var c = Math.Log(1 + nu);

            var thetaSquared = a * (b * c - 1);

            return thetaSquared; //in radians
        }
        public static double Vector3AngleBetween(Vector3 a, Vector3 b)
        {
            double dot = Vector3.Dot(Vector3.Normalize(a), Vector3.Normalize(b));
            return Math.Acos(Math.Clamp(dot, -1.0, 1.0));
        }

    }
}
