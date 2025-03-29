using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MCassignment_PHYS539_PeidusD_2025.Helpers
{
    public static class DoubleToFloatHelper
    {
        public static float[,,] ConvertToFloatArray(double[,,] input)
        {
            int d0 = input.GetLength(0);
            int d1 = input.GetLength(1);
            int d2 = input.GetLength(2);

            float[,,] output = new float[d0, d1, d2];

            for (int i = 0; i < d0; i++)
            {
                for (int j = 0; j < d1; j++)
                {
                    for (int k = 0; k < d2; k++)
                    {
                        output[i, j, k] = (float)input[i, j, k];
                    }
                }
            }
            return output;
        }
    }
    
}
