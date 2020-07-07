using System;
using System.Collections.Generic;

using ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;
using ENPC.Numerics;


namespace ENPC.NMontagne.Core.CoreFunctions.ChebyshevNets
{
    /// <summary>
    /// Class containing methods to generate Chebyshev nets on the unity sphere S².
    /// </summary>
    public static class ChebyshevOnUnitSphere
    {
        #region Primals Conditions

        /// <summary>
        /// Computes a Chebyshev net on the unity sphere from two primal conditions.
        /// </summary>
        /// <param name="normals_U"> The list of normals in the u-direction.</param>
        /// <param name="normals_V"> The list of normals in the v-direction.</param>
        /// <param name="mesh"> The mesh representing a Chebyshev net.</param>
        /// <exception cref="ArgumentNullException">  One of the list of points is empty.</exception>
        /// <exception cref="ArgumentException"> The first points of the lists must coincide.</exception>
        public static void Core_FromTwoPrimal(List<Point> normals_U, List<Point> normals_V, out HeMesh<Point> mesh)
        {
            #region Verifications

            if (normals_U.Count == 0) { throw new ArgumentNullException("The list of points in u-direction is empty."); }
            if (normals_V.Count == 0) { throw new ArgumentNullException("The list of points in v-direction is empty."); }
            if (normals_U[0].DistanceTo(normals_V[0]) > Settings._absolutePrecision)
            {
                throw new ArgumentException("The first points of the lists must coincide.");
            }

            #endregion

            #region Initilization

            int nb_U = normals_U.Count;
            int nb_V = normals_V.Count;

            mesh = new HeMesh<Point>();

            for (int i_U = 0; i_U < nb_U; i_U++)
            {
                mesh.AddVertex(normals_U[i_U]);
            }

            #endregion

            #region Propagation

            for (int i_V = 1; i_V < nb_V; i_V++)
            {
                mesh.AddVertex(normals_V[i_V]);
                for (int i_U = 1; i_U < nb_U; i_U++)
                {
                    // Get rotation axis
                    HeVertex<Point> v1 = mesh.GetVertex(i_U + ((i_V - 1) * nb_U));
                    HeVertex<Point> v2 = mesh.GetVertex((i_U - 1) + ((i_V) * nb_U));
                    Vector axis = (Vector)(v1.Position + v2.Position);
                    axis.Unitize();
                    // Get rotation quaternion
                    Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Get point to rotate
                    HeVertex<Point> v = mesh.GetVertex((i_U - 1) + ((i_V - 1) * nb_U));
                    Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                    // Compute the rotation
                    Quaternion qResult = qAxis * q * qAxis_1;

                    // Store result
                    HeVertex<Point> v12 = mesh.AddVertex(new Point(qResult.I, qResult.J, qResult.K));
                    mesh.AddFace(v, v1, v12, v2);
                }
            }

            #endregion
        }

        #endregion

        #region Dual Conditions

        /// <summary>
        /// Computes a Chebyshev net on the unity sphere from two open dual conditions.
        /// </summary>
        /// <param name="main"> The list of normals of the main diagonal.</param>
        /// <param name="minor"> The list of normals of a minor diagonal.</param>
        /// <param name="mesh"> The mesh representing a Chebyshev net.</param>
        /// <exception cref="ArgumentNullException"> One of the list of diagonal points is empty.</exception>
        /// <exception cref="ArgumentException"> For an open dual condition, the main diagonal must have one more element than minor diagonal.</exception>
        public static void Core_FromTwoOpenDual(List<Point> main, List<Point> minor, out HeMesh<Point> mesh)
        {
            #region Verifications

            if (main.Count == 0) { throw new ArgumentNullException("The list of points of the main diagonal is empty."); }
            if (minor.Count == 0) { throw new ArgumentNullException("The list of points of the main diagonal is empty."); }

            if (main.Count - 1 != minor.Count)
            {
                throw new ArgumentException("For an open dual condition, the main diagonal must have one more element than minor diagonal.");
            }

            #endregion

            #region Initialisation

            /* The main diagonal has i_Diag = 0 and the minor diagonal has i_Diag = 1
             * Thus, the diagonal following the minor diagonal have positive indices
             * while the diagonal preceding the main diagonal have negative indices.
             * 
             * This approach is not necesseraly convenient for the implementation but is more easily understood.
             * i_Diag is the index (positive or negative) of the diagonal, and i_Tab its corresponding index in the table
             * The rank of an element is its index in the diagonal.
             */

            int nb_Main = main.Count;

            int iMax_Diag = nb_Main - 1;
            int iMin_Diag = -iMax_Diag;

            // Local function to move from i_Diag to i_Tab
            int DiagtoTab(int i_Diag) => i_Diag + iMax_Diag;

            // Local function to get the number of element in a diagonal
            int RankCount(int i_Diag) => i_Diag < 0 ? nb_Main + i_Diag : nb_Main - i_Diag;

            /********** Initialization of the table **********/

            // Table to convert the pair (i_Diag, i_Rank) to the index of the vertex in the mesh
            int[][] converter = new int[2 * nb_Main - 1][];                                                 /**/
            for (int i_Diag = iMin_Diag; i_Diag < iMax_Diag + 1; i_Diag++)
            {
                converter[DiagtoTab(i_Diag)] = new int[RankCount(i_Diag)];
            }

            /********** Initialization of the mesh **********/

            mesh = new HeMesh<Point>();
            for (int i_Rank = 0; i_Rank < RankCount(0); i_Rank++)
            {
                mesh.AddVertex(main[i_Rank]);
                converter[DiagtoTab(0)][i_Rank] = i_Rank;
            }
            for (int i_Rank = 0; i_Rank < RankCount(1); i_Rank++)
            {
                mesh.AddVertex(minor[i_Rank]);
                converter[DiagtoTab(1)][i_Rank] = nb_Main + i_Rank;
            }

            #endregion

            #region Central Diagonal

            for (int i_Rank = 0; i_Rank < RankCount(-1); i_Rank++)
            {
                // Get rotation axis
                HeVertex<Point> v1 = mesh.GetVertex(converter[DiagtoTab(0)][i_Rank]);
                HeVertex<Point> v2 = mesh.GetVertex(converter[DiagtoTab(0)][i_Rank + 1]);
                Vector axis = (Vector)(v1.Position + v2.Position);
                axis.Unitize();
                // Get rotation quaternion
                Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                Quaternion qAxis_1 = qAxis.Inverse();
                // Get point to rotate
                HeVertex<Point> v = mesh.GetVertex(converter[DiagtoTab(1)][i_Rank]);
                Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                // Compute the rotation
                Quaternion qResult = qAxis * q * qAxis_1;

                // Store result
                HeVertex<Point> v12 = mesh.AddVertex(new Point(qResult.I, qResult.J, qResult.K));
                converter[DiagtoTab(-1)][i_Rank] = v12.Index;
                mesh.AddFace(v1, v12, v2, v);
            }

            #endregion

            #region Forward Propgation

            for (int i_Diag = 2; i_Diag < iMax_Diag + 1; i_Diag++)
            {
                for (int i_Rank = 0; i_Rank < RankCount(i_Diag); i_Rank++)
                {
                    // Get rotation axis
                    HeVertex<Point> v1 = mesh.GetVertex(converter[DiagtoTab(i_Diag - 1)][i_Rank]);
                    HeVertex<Point> v2 = mesh.GetVertex(converter[DiagtoTab(i_Diag - 1)][i_Rank + 1]);
                    Vector axis = (Vector)(v1.Position + v2.Position);
                    axis.Unitize();
                    // Get rotation quaternion
                    Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Get point to rotate
                    HeVertex<Point> v = mesh.GetVertex(converter[DiagtoTab(i_Diag - 2)][i_Rank + 1]);
                    Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                    // Compute the rotation
                    Quaternion qResult = qAxis * q * qAxis_1;

                    // Store result
                    HeVertex<Point> v12 = mesh.AddVertex(new Point(qResult.I, qResult.J, qResult.K));
                    converter[DiagtoTab(i_Diag)][i_Rank] = v12.Index;
                    mesh.AddFace(v1, v, v2, v12);
                }
            }

            #endregion

            #region Backward Propgation

            for (int i_Diag = -2; i_Diag > iMin_Diag - 1; i_Diag--)
            {
                for (int i_Rank = 0; i_Rank < RankCount(i_Diag); i_Rank++)
                {
                    // Get rotation axis
                    HeVertex<Point> v1 = mesh.GetVertex(converter[DiagtoTab(i_Diag + 1)][i_Rank]);
                    HeVertex<Point> v2 = mesh.GetVertex(converter[DiagtoTab(i_Diag + 1)][i_Rank + 1]);
                    Vector axis = (Vector)(v1.Position + v2.Position);
                    axis.Unitize();
                    // Get rotation quaternion
                    Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Get point to rotate
                    HeVertex<Point> v = mesh.GetVertex(converter[DiagtoTab(i_Diag + 2)][i_Rank + 1]);
                    Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                    // Compute the rotation
                    Quaternion qResult = qAxis * q * qAxis_1;

                    // Store result
                    HeVertex<Point> v12 = mesh.AddVertex(new Point(qResult.I, qResult.J, qResult.K));
                    converter[DiagtoTab(i_Diag)][i_Rank] = v12.Index;
                    mesh.AddFace(v1, v12, v2, v);
                }
            }

            #endregion

        }

        /// <summary>
        /// Computes a Chebyshev net on the unity sphere from two close dual conditions.
        /// </summary>
        /// <param name="main"> The list of normals of the main diagonal.</param>
        /// <param name="minor"> The list of normals of a minor diagonal.</param>
        /// <param name="iteration"> The number of iteration to compute.</param>
        /// <param name="mesh"> The mesh representing a Chebyshev net.</param>
        /// <exception cref="ArgumentNullException"> One of the list of diagonal points is empty.</exception>
        /// <exception cref="ArgumentException"> For an close dual condition, the main diagonal must have the same number of elements than minor diagonal.</exception>
        public static void Core_FromTwoCloseDual(List<Point> main, List<Point> minor, int iteration, out HeMesh<Point> mesh)
        {
            #region Verifications

            if (main.Count == 0) { throw new ArgumentNullException("The list of points of the main diagonal is empty."); }
            if (minor.Count == 0) { throw new ArgumentNullException("The list of points of the main diagonal is empty."); }

            if (main.Count != minor.Count)
            {
                throw new ArgumentException("For an close dual condition, the main diagonal must have the same number of elements than minor diagonal.");
            }

            #endregion

            #region Initialisation

            /* The main diagonal has i_Diag = 0 and the minor diagonal has i_Diag = 1
             * Thus, the diagonal following the minor diagonal have positive indices
             * while the diagonal preceding the main diagonal have negative indices.
             * 
             * This approach is not necesseraly convenient for the implementation but is more easily understood.
             * i_Diag is the index (positive or negative) of the diagonal, and i_Tab its corresponding index in the table
             * The rank of an element is its index in the diagonal.
             */

            int RankCount = main.Count;

            int iMax_Diag = iteration + 1;
            int iMin_Diag = -iMax_Diag;

            // Local function to move from i_Diag to i_Tab
            int DiagtoTab(int i_Diag) => i_Diag + iMax_Diag;

            // Local function to manage the looping rank of element in a diagonal
            int RankToTab(int i_Rank) => i_Rank == RankCount ? 0 : i_Rank;

            /********** Initialization of the table **********/

            // Table to convert the pair (i_Diag, i_Rank) to the index of the vertex in the mesh
            int[,] converter = new int[2 * iMax_Diag + 1, RankCount];

            /********** Initialization of the mesh **********/

            mesh = new HeMesh<Point>();
            for (int i_Rank = 0; i_Rank < RankCount; i_Rank++)
            {
                mesh.AddVertex(main[i_Rank]);
                converter[DiagtoTab(0), i_Rank] = i_Rank;
            }
            for (int i_Rank = 0; i_Rank < RankCount; i_Rank++)
            {
                mesh.AddVertex(minor[i_Rank]);
                converter[DiagtoTab(1), i_Rank] = RankCount + i_Rank;
            }

            #endregion

            #region Central Diagonal

            for (int i_Rank = 0; i_Rank < RankCount; i_Rank++)
            {
                // Get rotation axis
                HeVertex<Point> v1 = mesh.GetVertex(converter[DiagtoTab(0),i_Rank]);
                HeVertex<Point> v2 = mesh.GetVertex(converter[DiagtoTab(0), RankToTab(i_Rank + 1)]);
                Vector axis = (Vector)(v1.Position + v2.Position);
                axis.Unitize();
                // Get rotation quaternion
                Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                Quaternion qAxis_1 = qAxis.Inverse();
                // Get point to rotate
                HeVertex<Point> v = mesh.GetVertex(converter[DiagtoTab(1),i_Rank]);
                Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                // Compute the rotation
                Quaternion qResult = qAxis * q * qAxis_1;

                // Store result
                HeVertex<Point> v12 = mesh.AddVertex(new Point(qResult.I, qResult.J, qResult.K));
                converter[DiagtoTab(-1),i_Rank] = v12.Index;
                mesh.AddFace(v1, v12, v2, v);
            }

            #endregion

            #region Forward Propgation

            for (int i_Diag = 2; i_Diag < iMax_Diag + 1; i_Diag++)
            {
                for (int i_Rank = 0; i_Rank < RankCount; i_Rank++)
                {
                    // Get rotation axis
                    HeVertex<Point> v1 = mesh.GetVertex(converter[DiagtoTab(i_Diag - 1), i_Rank]);
                    HeVertex<Point> v2 = mesh.GetVertex(converter[DiagtoTab(i_Diag - 1), RankToTab(i_Rank + 1)]);
                    Vector axis = (Vector)(v1.Position + v2.Position);
                    axis.Unitize();
                    // Get rotation quaternion
                    Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Get point to rotate
                    HeVertex<Point> v = mesh.GetVertex(converter[DiagtoTab(i_Diag - 2), RankToTab(i_Rank + 1)]);
                    Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                    // Compute the rotation
                    Quaternion qResult = qAxis * q * qAxis_1;

                    // Store result
                    HeVertex<Point> v12 = mesh.AddVertex(new Point(qResult.I, qResult.J, qResult.K));
                    converter[DiagtoTab(i_Diag), i_Rank] = v12.Index;
                    mesh.AddFace(v1, v, v2, v12);
                }
            }

            #endregion

            #region Backward Propgation

            for (int i_Diag = -2; i_Diag > iMin_Diag - 1; i_Diag--)
            {
                for (int i_Rank = 0; i_Rank < RankCount; i_Rank++)
                {
                    // Get rotation axis
                    HeVertex<Point> v1 = mesh.GetVertex(converter[DiagtoTab(i_Diag + 1), i_Rank]);
                    HeVertex<Point> v2 = mesh.GetVertex(converter[DiagtoTab(i_Diag + 1), RankToTab(i_Rank + 1)]);
                    Vector axis = (Vector)(v1.Position + v2.Position);
                    axis.Unitize();
                    // Get rotation quaternion
                    Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Get point to rotate
                    HeVertex<Point> v = mesh.GetVertex(converter[DiagtoTab(i_Diag + 2), RankToTab(i_Rank + 1)]);
                    Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                    // Compute the rotation
                    Quaternion qResult = qAxis * q * qAxis_1;

                    // Store result
                    HeVertex<Point> v12 = mesh.AddVertex(new Point(qResult.I, qResult.J, qResult.K));
                    converter[DiagtoTab(i_Diag), i_Rank] = v12.Index;
                    mesh.AddFace(v1, v12, v2, v);
                }
            }

            #endregion
        }

        #endregion

        #region Singularity Conditions

        /// <summary>
        /// Computes a Chebyshev net on the unity sphere from a rosette singularity.
        /// </summary>
        /// <param name="centre"> The singularity point.</param>
        /// <param name="ring"> The list of normals of the ring following the singularity.</param>
        /// <param name="iteration"> The number of iteration to compute.</param>
        /// <param name="mesh"> The mesh representing a Chebyshev net.</param>
        /// <exception cref="ArgumentNullException"> The list of points of the ring is empty.</exception>
        public static void Core_Rosette(Point centre, List<Point> ring, int iteration, out HeMesh<Point> mesh)
        {
            #region Verifications

            if (ring.Count == 0) { throw new ArgumentNullException("The list of points of the ring is empty."); }

            #endregion

            #region Initialisation

            /* The main diagonal has i_Diag = 0 and the minor diagonal has i_Diag = 1
             * Thus, the diagonal following the minor diagonal have positive indices
             * while the diagonal preceding the main diagonal have negative indices.
             * 
             * This approach is not necesseraly convenient for the implementation but is more easily understood.
             * i_Diag is the index (positive or negative) of the diagonal, and i_Tab its corresponding index in the table
             * The rank of an element is its index in the diagonal.
             */

            int RankCount = ring.Count;

            int iMax_Diag = iteration + 1;

            // Local function to manage the looping rank of element in a diagonal
            int RankToTab(int i_Rank) => i_Rank == RankCount ? 0 : i_Rank;

            /********** Initialization of the table **********/

            // Table to convert the pair (i_Diag, i_Rank) to the index of the vertex in the mesh
            int[,] converter = new int[iMax_Diag + 1, RankCount];

            /********** Initialization of the mesh **********/

            mesh = new HeMesh<Point>();

            mesh.AddVertex(centre);

            for (int i_Rank = 0; i_Rank < RankCount; i_Rank++)
            {
                converter[0, i_Rank] = 0;

                mesh.AddVertex(ring[i_Rank]);
                converter[1, i_Rank] = 1 + i_Rank;
            }

            #endregion

            #region Propgation

            for (int i_Diag = 2; i_Diag < iMax_Diag + 1; i_Diag++)
            {
                for (int i_Rank = 0; i_Rank < RankCount; i_Rank++)
                {
                    // Get rotation axis
                    HeVertex<Point> v1 = mesh.GetVertex(converter[i_Diag - 1, i_Rank]);
                    HeVertex<Point> v2 = mesh.GetVertex(converter[i_Diag - 1, RankToTab(i_Rank + 1)]);
                    Vector axis = (Vector)(v1.Position + v2.Position);
                    axis.Unitize();
                    // Get rotation quaternion
                    Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Get point to rotate
                    HeVertex<Point> v = mesh.GetVertex(converter[i_Diag - 2, RankToTab(i_Rank + 1)]);
                    Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                    // Compute the rotation
                    Quaternion qResult = qAxis * q * qAxis_1;

                    // Store result
                    HeVertex<Point> v12 = mesh.AddVertex(new Point(qResult.I, qResult.J, qResult.K));
                    converter[i_Diag, i_Rank] = v12.Index;
                    mesh.AddFace(v1, v, v2, v12);
                }
            }

            #endregion
        }

        #endregion

    }
}
