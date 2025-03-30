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

        // Prepare a buffer for thread-safe merging
        double[,,] accumulatedResult = new double[nx, ny, nz];

        Parallel.For(0, nx, ix =>
        {
            double[,,] localBuffer = new double[nx, ny, nz];

            for (int iy = 0; iy < ny; iy++)
            {
                for (int iz = 0; iz < nz; iz++)
                {
                    double value = photonMapDoseMatrix[ix, iy, iz];
                    if (value == 0) continue;

                    for (int dx = 0; dx < kx; dx++)
                    {
                        int tx = ix + dx - cx;
                        if (tx < 0 || tx >= nx) continue;

                        for (int dy = 0; dy < ky; dy++)
                        {
                            int ty = iy + dy - cy;
                            if (ty < 0 || ty >= ny) continue;

                            for (int dz = 0; dz < kz; dz++)
                            {
                                int tz = iz + dz - cz;
                                if (tz < 0 || tz >= nz) continue;

                                double dose = value * kernelDoseMatrix[dx, dy, dz];
                                if (dose == 0) continue;

                                localBuffer[tx, ty, tz] += dose;
                            }
                        }
                    }
                }
            }

            lock (accumulatedResult)
            {
                for (int x = 0; x < nx; x++)
                    for (int y = 0; y < ny; y++)
                        for (int z = 0; z < nz; z++)
                            photonMap.DepositEnergy(x, y, z, localBuffer[x, y, z]) ;
            }
        });

        double maxDose = photonMapDoseMatrix.Cast<double>().Max();
        ExportPeakSlices(photonMapDoseMatrix, "TEST");
        Console.WriteLine($"Brute-force convolution completed. Max dose: {maxDose:G6}");
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
