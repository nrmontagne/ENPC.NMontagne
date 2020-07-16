using System;
using System.Collections.Generic;

using Euc = ENPC.Geometry.Euclidean;
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
        public static void Core_FromTwoPrimal(List<Euc.Point> normals_U, List<Euc.Point> normals_V, out HeMesh<Euc.Point> mesh)
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

            mesh = new HeMesh<Euc.Point>();

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
                    HeVertex<Euc.Point> v1 = mesh.GetVertex(i_U + ((i_V - 1) * nb_U));
                    HeVertex<Euc.Point> v2 = mesh.GetVertex((i_U - 1) + ((i_V) * nb_U));
                    Euc.Vector axis = (Euc.Vector)(v1.Position + v2.Position);
                    axis.Unitize();
                    // Get rotation quaternion
                    Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Get point to rotate
                    HeVertex<Euc.Point> v = mesh.GetVertex((i_U - 1) + ((i_V - 1) * nb_U));
                    Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                    // Compute the rotation
                    Quaternion qResult = qAxis * q * qAxis_1;

                    // Store result
                    HeVertex<Euc.Point> v12 = mesh.AddVertex(new Euc.Point(qResult.I, qResult.J, qResult.K));
                    mesh.AddFace(v, v1, v12, v2);
                }
            }

            #endregion
        }

        /// <summary>
        /// Computes two patched Chebyshev nets on the unity sphere from three primal conditions.
        /// </summary>
        /// <param name="normals_A"> The list of normals for the first Chebyshev net.</param>
        /// <param name="normals_B"> The list of normals common to the two Chebyshev nets.</param>
        /// <param name="normals_C"> The list of normals for the second Chebyshev net.</param>
        /// <param name="mesh"> The mesh representing the two patched Chebyshev nets.</param>
        /// <exception cref="ArgumentNullException">  One of the list of points is empty.</exception>
        /// <exception cref="ArgumentException"> The first points of the lists must coincide.</exception>
        public static void Core_FromThreePrimal(List<Euc.Point> normals_A, List<Euc.Point> normals_B, List<Euc.Point> normals_C, out HeMesh<Euc.Point> mesh)
        {
            #region Verifications

            if (normals_A.Count == 0) { throw new ArgumentNullException("The list of points in A is empty."); }
            if (normals_B.Count == 0) { throw new ArgumentNullException("The list of points in B is empty."); }
            if (normals_C.Count == 0) { throw new ArgumentNullException("The list of points in C is empty."); }
            if (normals_A[0].DistanceTo(normals_B[0]) > Settings._absolutePrecision || normals_A[0].DistanceTo(normals_C[0]) > Settings._absolutePrecision)
            {
                throw new ArgumentException("The first points of the lists must coincide.");
            }

            #endregion

            #region Initilization

            int nb_A = normals_A.Count;
            int nb_B = normals_B.Count;
            int nb_C = normals_C.Count;

            mesh = new HeMesh<Euc.Point>();

            for (int i_U = 0; i_U < nb_B; i_U++)
            {
                mesh.AddVertex(normals_B[i_U]);
            }

            #endregion

            #region Propagation BA

            for (int i_V = 1; i_V < nb_A; i_V++)
            {
                mesh.AddVertex(normals_A[i_V]);
                for (int i_U = 1; i_U < nb_B; i_U++)
                {
                    // Get rotation axis
                    HeVertex<Euc.Point> v1 = mesh.GetVertex(i_U + ((i_V - 1) * nb_B));
                    HeVertex<Euc.Point> v2 = mesh.GetVertex((i_U - 1) + ((i_V) * nb_B));
                    Euc.Vector axis = (Euc.Vector)(v1.Position + v2.Position);
                    axis.Unitize();
                    // Get rotation quaternion
                    Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Get point to rotate
                    HeVertex<Euc.Point> v = mesh.GetVertex((i_U - 1) + ((i_V - 1) * nb_B));
                    Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                    // Compute the rotation
                    Quaternion qResult = qAxis * q * qAxis_1;

                    // Store result
                    HeVertex<Euc.Point> v12 = mesh.AddVertex(new Euc.Point(qResult.I, qResult.J, qResult.K));
                    mesh.AddFace(v, v1, v12, v2);
                }
            }

            #endregion

            #region Propagation BC

            int nb_VertexAB = nb_A * nb_B;
            for (int i_V = 1; i_V < nb_A; i_V++)
            {
                mesh.AddVertex(normals_A[i_V]);
                for (int i_U = 1; i_U < nb_B; i_U++)
                {
                    // Manage the vertex indices
                    int v_Index = i_V == 1 ? i_U - 1 : nb_VertexAB + (i_U - 1) + ((i_V - 1) * nb_B);
                    int v1_Index = i_V == 1 ? i_U : nb_VertexAB + i_U + ((i_V - 1) * nb_B);
                    int v2_Index = nb_VertexAB + (i_U - 1) + ((i_V) * nb_B);
                    int v12_Index = nb_VertexAB + i_U + ((i_V) * nb_B);

                    // Get rotation axis
                    HeVertex<Euc.Point> v1 = mesh.GetVertex(v1_Index);
                    HeVertex<Euc.Point> v2 = mesh.GetVertex(v2_Index);
                    Euc.Vector axis = (Euc.Vector)(v1.Position + v2.Position);
                    axis.Unitize();
                    // Get rotation quaternion
                    Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Get point to rotate
                    HeVertex<Euc.Point> v = mesh.GetVertex(v_Index);
                    Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                    // Compute the rotation
                    Quaternion qResult = qAxis * q * qAxis_1;

                    // Store result
                    HeVertex<Euc.Point> v12 = mesh.AddVertex(new Euc.Point(qResult.I, qResult.J, qResult.K));
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
        public static void Core_FromTwoOpenDual(List<Euc.Point> main, List<Euc.Point> minor, out HeMesh<Euc.Point> mesh)
        {
            #region Verifications

            if (main.Count == 0) { throw new ArgumentNullException("The list of points of the main diagonal is empty."); }
            if (minor.Count == 0) { throw new ArgumentNullException("The list of points of the minor diagonal is empty."); }

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

            mesh = new HeMesh<Euc.Point>();
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
                HeVertex<Euc.Point> v1 = mesh.GetVertex(converter[DiagtoTab(0)][i_Rank]);
                HeVertex<Euc.Point> v2 = mesh.GetVertex(converter[DiagtoTab(0)][i_Rank + 1]);
                Euc.Vector axis = (Euc.Vector)(v1.Position + v2.Position);
                axis.Unitize();
                // Get rotation quaternion
                Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                Quaternion qAxis_1 = qAxis.Inverse();
                // Get point to rotate
                HeVertex<Euc.Point> v = mesh.GetVertex(converter[DiagtoTab(1)][i_Rank]);
                Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                // Compute the rotation
                Quaternion qResult = qAxis * q * qAxis_1;

                // Store result
                HeVertex<Euc.Point> v12 = mesh.AddVertex(new Euc.Point(qResult.I, qResult.J, qResult.K));
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
                    HeVertex<Euc.Point> v1 = mesh.GetVertex(converter[DiagtoTab(i_Diag - 1)][i_Rank]);
                    HeVertex<Euc.Point> v2 = mesh.GetVertex(converter[DiagtoTab(i_Diag - 1)][i_Rank + 1]);
                    Euc.Vector axis = (Euc.Vector)(v1.Position + v2.Position);
                    axis.Unitize();
                    // Get rotation quaternion
                    Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Get point to rotate
                    HeVertex<Euc.Point> v = mesh.GetVertex(converter[DiagtoTab(i_Diag - 2)][i_Rank + 1]);
                    Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                    // Compute the rotation
                    Quaternion qResult = qAxis * q * qAxis_1;

                    // Store result
                    HeVertex<Euc.Point> v12 = mesh.AddVertex(new Euc.Point(qResult.I, qResult.J, qResult.K));
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
                    HeVertex<Euc.Point> v1 = mesh.GetVertex(converter[DiagtoTab(i_Diag + 1)][i_Rank]);
                    HeVertex<Euc.Point> v2 = mesh.GetVertex(converter[DiagtoTab(i_Diag + 1)][i_Rank + 1]);
                    Euc.Vector axis = (Euc.Vector)(v1.Position + v2.Position);
                    axis.Unitize();
                    // Get rotation quaternion
                    Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Get point to rotate
                    HeVertex<Euc.Point> v = mesh.GetVertex(converter[DiagtoTab(i_Diag + 2)][i_Rank + 1]);
                    Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                    // Compute the rotation
                    Quaternion qResult = qAxis * q * qAxis_1;

                    // Store result
                    HeVertex<Euc.Point> v12 = mesh.AddVertex(new Euc.Point(qResult.I, qResult.J, qResult.K));
                    converter[DiagtoTab(i_Diag)][i_Rank] = v12.Index;
                    mesh.AddFace(v1, v12, v2, v);
                }
            }

            #endregion

        }

        /// <summary>
        /// Computes two patched Chebyshev net on the unity sphere from three dual conditions.
        /// </summary>
        /// <param name="main_A"> The list of normals of the main diagonal for the first Chebyshev.</param>
        /// <param name="minor_A"> The list of normals of a minor diagonal for the first Chebyshev.</param>
        /// <param name="main_B">  The list of normals of the main diagonal for the second Chebyshev.</param>
        /// <param name="mesh"> The mesh representing a the two patched Chebyshev nets.</param>
        [Obsolete]
        public static void Core_FromThreeOpenDual(List<Euc.Point> main_A, List<Euc.Point> minor_A, List<Euc.Point> main_B, out HeMesh<Euc.Point> mesh)
        {

            #region Verifications

            if (main_A is null || minor_A is null || main_B is null) { throw new ArgumentNullException(); }
            if (main_A.Count == 0 || minor_A.Count == 0 || main_B.Count == 0) { throw new ArgumentException("The list of normals shouldn't be empty"); }
            if (!main_A[0].Equals(main_B[0])) { throw new ArgumentException("The first normals in Normals_A_0 and Normals_B_0 should coincide."); }
            if (main_A.Count - 1 != minor_A.Count) { throw new ArgumentException("Normals_A_0 must have one more element than Normals_A_1."); }
            if (main_B.Count > main_A.Count)
            {
                throw new ArgumentException("For now, the number of elements in Normals_B_0 should exceed the number of elements in Normals_A_0.");
            }

            #endregion

            #region Utilities

            int rowCount = main_A.Count + main_B.Count - 1;

            int RankCount(int row) => row < main_A.Count ? row + 1 : main_A.Count - (row - main_A.Count + 1);

            #endregion

            #region Initialization

            Euc.Point[][] normals = new Euc.Point[rowCount][];
            for (int row = 0; row < rowCount; row++)
            {
                normals[row] = new Euc.Point[RankCount(row)];

                if (row < main_A.Count)
                {
                    normals[row][0] = main_A[main_A.Count - 1 - row];
                    if (row > 0) { normals[row][1] = minor_A[minor_A.Count - row]; }
                }
                else { normals[row][0] = main_B[row - (main_A.Count - 1)]; }
            }

            #endregion

            #region Propagation Forward

            /******************** Propagation for Normals A ********************/

            for (int row = 2; row < main_A.Count; row++)
            {
                for (int rank = 2; rank < RankCount(row); rank++)
                {
                    // Rotation Quaternions
                    Euc.Vector axis = (Euc.Vector) (normals[row][rank - 1] + normals[row - 1][rank - 1]);
                    axis.Unitize();
                    double angle = Math.PI;
                    Quaternion qAxis = Quaternion.Unit(axis, angle);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Vector Quaternion
                    Quaternion qVector = new Quaternion(0, normals[row - 1][rank - 2].X, normals[row - 1][rank - 2].Y, normals[row - 1][rank - 2].Z);
                    // Compute rotation about "axis" of "angle" of "vector"
                    Quaternion qResult = qAxis * qVector * qAxis_1;
                    normals[row][rank] = new Euc.Point(qResult.I, qResult.J, qResult.K);
                }
            }

            /******************** Propagation for Normals B ********************/

            for (int row = main_A.Count; row < rowCount; row++)
            {
                for (int rank = 1; rank < RankCount(row); rank++)
                {
                    // Rotation Quaternions
                    Euc.Vector axis = (Euc.Vector) (normals[row][rank - 1] + normals[row - 1][rank + 1]);
                    axis.Unitize();
                    double angle = Math.PI;
                    Quaternion qAxis = Quaternion.Unit(axis, angle);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Vector Quaternion
                    Quaternion qVector = new Quaternion(0, normals[row - 1][rank].X, normals[row - 1][rank].Y, normals[row - 1][rank].Z);
                    // Compute rotation about "axis" of "angle" of "vector"
                    Quaternion qResult = qAxis * qVector * qAxis_1;
                    normals[row][rank] = new Euc.Point(qResult.I, qResult.J, qResult.K);
                }
            }

            #endregion

            #region Generate HeMesh

            mesh = new HeMesh<Euc.Point>();

            // Create a converter that will convert the indices of a normal in the array of vertices
            // (previously computed) to the index of the same normal in the list of vertices of the new mesh.

            int[][] conv = new int[rowCount][];
            for (int row = 0; row < rowCount; row++)
            {
                conv[row] = new int[RankCount(row)];
            }

            /******************** Add Vertices ********************/

            for (int row = 0; row < rowCount; row++)
            {
                for (int rank = 0; rank < RankCount(row); rank++)
                {
                    conv[row][rank] = mesh.AddVertex((Euc.Point)normals[row][rank]).Index;
                }
            }

            /******************** Add Quad Faces from Normals A ********************/

            for (int row = 1; row < main_A.Count - 1; row++)
            {
                for (int rank = 0; rank < RankCount(row) - 1; rank++)
                {
                    mesh.AddFace(mesh.GetVertex(conv[row][rank]), mesh.GetVertex(conv[row][rank + 1]),
                        mesh.GetVertex(conv[row + 1][rank + 2]), mesh.GetVertex(conv[row + 1][rank + 1]));
                }
            }

            /******************** Add Quad Faces from Normals B ********************/

            for (int row = main_A.Count - 1; row < rowCount - 2; row++)
            {
                for (int rank = 1; rank < RankCount(row) - 1; rank++)
                {
                    mesh.AddFace(mesh.GetVertex(conv[row][rank]), mesh.GetVertex(conv[row][rank + 1]),
                        mesh.GetVertex(conv[row + 1][rank]), mesh.GetVertex(conv[row + 1][rank - 1]));
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
        public static void Core_FromTwoCloseDual(List<Euc.Point> main, List<Euc.Point> minor, int iteration, out HeMesh<Euc.Point> mesh)
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

            mesh = new HeMesh<Euc.Point>();
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
                HeVertex<Euc.Point> v1 = mesh.GetVertex(converter[DiagtoTab(0),i_Rank]);
                HeVertex<Euc.Point> v2 = mesh.GetVertex(converter[DiagtoTab(0), RankToTab(i_Rank + 1)]);
                Euc.Vector axis = (Euc.Vector)(v1.Position + v2.Position);
                axis.Unitize();
                // Get rotation quaternion
                Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                Quaternion qAxis_1 = qAxis.Inverse();
                // Get point to rotate
                HeVertex<Euc.Point> v = mesh.GetVertex(converter[DiagtoTab(1),i_Rank]);
                Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                // Compute the rotation
                Quaternion qResult = qAxis * q * qAxis_1;

                // Store result
                HeVertex<Euc.Point> v12 = mesh.AddVertex(new Euc.Point(qResult.I, qResult.J, qResult.K));
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
                    HeVertex<Euc.Point> v1 = mesh.GetVertex(converter[DiagtoTab(i_Diag - 1), i_Rank]);
                    HeVertex<Euc.Point> v2 = mesh.GetVertex(converter[DiagtoTab(i_Diag - 1), RankToTab(i_Rank + 1)]);
                    Euc.Vector axis = (Euc.Vector)(v1.Position + v2.Position);
                    axis.Unitize();
                    // Get rotation quaternion
                    Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Get point to rotate
                    HeVertex<Euc.Point> v = mesh.GetVertex(converter[DiagtoTab(i_Diag - 2), RankToTab(i_Rank + 1)]);
                    Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                    // Compute the rotation
                    Quaternion qResult = qAxis * q * qAxis_1;

                    // Store result
                    HeVertex<Euc.Point> v12 = mesh.AddVertex(new Euc.Point(qResult.I, qResult.J, qResult.K));
                    converter[DiagtoTab(i_Diag), i_Rank] = v12.Index;
                    mesh.AddFace(v1, v, v2, v12);
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
        public static void Core_Rosette(Euc.Point centre, List<Euc.Point> ring, int iteration, out HeMesh<Euc.Point> mesh)
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

            mesh = new HeMesh<Euc.Point>();

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
                    HeVertex<Euc.Point> v1 = mesh.GetVertex(converter[i_Diag - 1, i_Rank]);
                    HeVertex<Euc.Point> v2 = mesh.GetVertex(converter[i_Diag - 1, RankToTab(i_Rank + 1)]);
                    Euc.Vector axis = (Euc.Vector)(v1.Position + v2.Position);
                    axis.Unitize();
                    // Get rotation quaternion
                    Quaternion qAxis = Quaternion.Unit(axis, Math.PI);
                    Quaternion qAxis_1 = qAxis.Inverse();
                    // Get point to rotate
                    HeVertex<Euc.Point> v = mesh.GetVertex(converter[i_Diag - 2, RankToTab(i_Rank + 1)]);
                    Quaternion q = new Quaternion(0, v.Position.X, v.Position.Y, v.Position.Z);
                    // Compute the rotation
                    Quaternion qResult = qAxis * q * qAxis_1;

                    // Store result
                    HeVertex<Euc.Point> v12 = mesh.AddVertex(new Euc.Point(qResult.I, qResult.J, qResult.K));
                    converter[i_Diag, i_Rank] = v12.Index;
                    mesh.AddFace(v1, v, v2, v12);
                }
            }

            #endregion
        }

        #endregion

    }
}
