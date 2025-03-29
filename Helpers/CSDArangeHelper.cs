using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MCassignment_PHYS539_PeidusD_2025.Helpers
{
    public class CSDArangeHelper
    {
        private static readonly (double E, double rangeCm)[] CSDArangeTable = new[]
        {
            (0.1,0.014),
            (0.5,0.18),
            (1.0,0.44),
            (2.0,1.0),
            (4.0,2.0),
            (10.0,5.0)
        };

        public static double GetCSDArange(double eMeV)
        {
            for (int i=0;i<CSDArangeTable.Length-1;i++)
            {
                var (E1,R1) = CSDArangeTable[i];
                var (E2,R2) = CSDArangeTable[i+1];

                if (eMeV >= E1&& eMeV <= E2)
                {
                    double fraction = (eMeV - E1) / (E2 - E1);
                    var stopPower = R1 + fraction * (R2 - R1);
                    return stopPower;
                }
            }

            if (eMeV < CSDArangeTable[0].E)
                return CSDArangeTable[0].rangeCm * (eMeV / CSDArangeTable[0].E);
            if (eMeV < CSDArangeTable[^1].E)
                return CSDArangeTable[^1].rangeCm * (eMeV / CSDArangeTable[^1].E);
            
            return 0.01; //fallback to prevent division by 0
        }
        public static double GetStoppingPower (double eMeV)
        {
            double range = GetCSDArange(eMeV);
            return eMeV / range;
        }
    }
}
