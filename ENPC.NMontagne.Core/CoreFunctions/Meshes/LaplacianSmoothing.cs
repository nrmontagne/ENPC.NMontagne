using System;
using System.Collections.Generic;

using Euc = ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;


namespace ENPC.NMontagne.Core.CoreFunctions.Meshes
{
    /// <summary>
    /// Class containing methods to compute the laplcian smoothing of a mesh.
    /// </summary>
    public static class LaplacianSmoothing
    {
        /// <summary>
        /// Performs the Laplacian smoothing of a mesh.
        /// </summary>
        /// <param name="mesh"> The mesh to operate on.</param>
        /// <param name="iteration"> The number of smoothing iterations.</param>
        /// <param name="condition"> The boundary condition : <br/> 0 : free edges; 1 : fixed boundary;</param>
        /// <param name="defMEsh"> The smoothed mesh.</param>
        public static void Core_NotWeighted(HeMesh<Euc.Point> mesh, int iteration, int condition, out HeMesh<Euc.Point> defMEsh)
        {
            defMEsh = (HeMesh<Euc.Point>)mesh.Clone();
            defMEsh.LaplacianSmoothing(0.5, iteration, condition);
        }
    }
}
