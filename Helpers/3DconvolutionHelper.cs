using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Particles;
using System;
using System.Linq;
using System.Threading.Tasks;

public static class ManualConvolver3D
{
    public static void Convolve(LargeWaterTank photonMap, WaterTankPhantom_UniformResolution kernel)
    {
        Console.WriteLine("Brute-force convolution started...");

        int nx = photonMap.Nx;
        int ny = photonMap.Ny;
        int nz = photonMap.Nz;

        int kx = kernel.Nx;
        int ky = kernel.Ny;
        int kz = kernel.Nz;

        var kernelDoseMatrix = kernel.GetDoseGrid();
        var photonMapDoseMatrix = photonMap.GetDoseGrid();

        // Find kernel center as the coordinates of the max value
        int cx = 0, cy = 0, cz = 0;
        double maxVal = double.MinValue;

        for (int i = 0; i < kx; i++)
            for (int j = 0; j < ky; j++)
                for (int k = 0; k < kz; k++)
                {
                    if (kernelDoseMatrix[i, j, k] > maxVal)
                    {
                        maxVal = kernelDoseMatrix[i, j, k];
                        cx = i;
                        cy = j;
                        cz = k;
                    }
                }
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                for (int iz = 0; iz < ny; iz++)
                {
                    for (int i0 = 0; i0 < kx; i0++)
                    {
                        for (int i1 = 0; i1 < ky; i1++)
                        {
                            for (int i2 = 0; i2 < kz; i2++)
                            {
                                if (i0 < kx && i1 < ky && i2 < kz)
                                {
                                    var test1 = photonMapDoseMatrix[ix, iy, iz];
                                    var test2 = kernelDoseMatrix[i0, i1, i2];
                                    var test3 = photonMapDoseMatrix[i0, i1, i2];
                                    if (test2 != 0)
                                    {
                                        photonMapDoseMatrix[i0, i1, i2] += test1 * test2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        

        
        double maxDose = photonMapDoseMatrix.Cast<double>().Max();
        Console.WriteLine($"Brute-force convolution completed. Max dose: {maxDose:G6}");
        ExportPeakSlices(photonMapDoseMatrix, "DoseMatrix");

    }
    public static void ExportPeakSlices(double[,,] doseGrid, string baseName = "dose")
    {
        int nx = doseGrid.GetLength(0);
        int ny = doseGrid.GetLength(1);
        int nz = doseGrid.GetLength(2);

        // --- Find Z index with max XY slice ---
        double maxSumZ = double.MinValue;
        int maxZ = 0;

        for (int z = 0; z < nz; z++)
        {
            double sum = 0;
            for (int x = 0; x < nx; x++)
                for (int y = 0; y < ny; y++)
                    sum += doseGrid[x, y, z];

            if (sum > maxSumZ)
            {
                maxSumZ = sum;
                maxZ = z;
            }
        }

        // --- Find Y index with max XZ slice ---
        double maxSumY = double.MinValue;
        int maxY = 0;

        for (int y = 0; y < ny; y++)
        {
            double sum = 0;
            for (int x = 0; x < nx; x++)
                for (int z = 0; z < nz; z++)
                    sum += doseGrid[x, y, z];

            if (sum > maxSumY)
            {
                maxSumY = sum;
                maxY = y;
            }
        }

        // --- Export XY @ maxZ ---
        using (var writer = new StreamWriter($"{baseName}_XY_z{maxZ}.csv"))
        {
            for (int y = 0; y < ny; y++)
            {
                for (int x = 0; x < nx; x++)
                {
                    writer.Write(doseGrid[x, y, maxZ].ToString("G6"));
                    if (x < nx - 1) writer.Write(";");
                }
                writer.WriteLine();
            }
        }

        // --- Export XZ @ maxY ---
        using (var writer = new StreamWriter($"{baseName}_XZ_y{maxY}.csv"))
        {
            for (int z = 0; z < nz; z++)
            {
                for (int x = 0; x < nx; x++)
                {
                    writer.Write(doseGrid[x, maxY, z].ToString("G6"));
                    if (x < nx - 1) writer.Write(";");
                }
                writer.WriteLine();
            }
        }

        Console.WriteLine($"✅ Saved XY slice @ Z={maxZ} and XZ slice @ Y={maxY}");
    }
}
