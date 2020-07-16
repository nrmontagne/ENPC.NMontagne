using System;
using System.Linq;
using System.Collections.Generic;

using Acc = Accord.Math;

using Euc = ENPC.Geometry.Euclidean;
using LinAlg = ENPC.Numerics.LinearAlgebra;

using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;


namespace ENPC.NMontagne.Core.CoreFunctions.VossNets
{
    /// <summary>
    /// Class containing methods to generate Voss nets from its Gauss map using linear optimisation.
    /// </summary>
    public static class GaussToVoss_LinOpt
    {
        /******************** General Voss Net ********************/

        /// <summary> 
        /// Computes Voss net from a Gauss map using a target length for the edges.
        /// </summary>
        /// <param name="gaussMap"> The Gauss map of the Voss net to generate.</param>
        /// <param name="Kp"> Evaluates whether the Gaussian curvature should be positive or not.</param>
        /// <param name="vossNet"> The Voss net generated from the input Gauss map.</param>
        public static void Core_TargetLength_Acc(HeMesh<Euc.Point> gaussMap, bool Kp, out HeMesh<Euc.Point> vossNet)
        {

            #region Variables Declaration

            vossNet = new HeMesh<Euc.Point>();

            int nbVertex;       // Equal to the number of faces in the Gauss image.
            int nbEdge;         // Equal to the number of internal edges in the Gauss image.
            int nbFace;         // Equal to the number of internal vertices in the Gauss image. 

            List<int> e_VtoG = new List<int>();     // Convert Voss edge indices into Gauss Halfedge indices
            List<int> f_VtoG = new List<int>();     // Convert Voss face indices into Gauss vertex indices

            #endregion

            #region Variables Initialisation

            // Number of Voss vertices
            nbVertex = gaussMap.FaceCount;

            // Number of Voss edges
            for (int i = 0; i < gaussMap.EdgeCount / 2; i++)
            {
                HeEdge<Euc.Point> edge = gaussMap.GetEdge(2 * i);
                if (edge.IsBoundary() || edge.PairEdge.IsBoundary()) { continue; }
                e_VtoG.Add(2 * i);
            }
            nbEdge = e_VtoG.Count;

            //Number of Voss faces
            for (int i = 0; i < gaussMap.VertexCount; i++)
            {
                if (gaussMap.GetVertex(i).IsBoundary()) { continue; }
                f_VtoG.Add(i);
            }
            nbFace = f_VtoG.Count;

            #endregion

            #region Edge Directions 

            Euc.Vector[] dir_Edges = new Euc.Vector[nbEdge];
            for (int i_Edge = 0; i_Edge < nbEdge; i_Edge++)
            {
                Euc.Vector dir = Euc.Vector.CrossProduct((Euc.Vector)gaussMap.GetEdge(e_VtoG[i_Edge]).EndVertex.Position,
                    (Euc.Vector)gaussMap.GetEdge(e_VtoG[i_Edge]).StartVertex.Position);

                dir_Edges[i_Edge] = dir;
            }

            #endregion


            #region Edges from Vertices : EfromV

            // With this matrix we have E = EfromV * V.

            double[,] EfromV = Acc.Matrix.Create(3 * nbEdge, 3 * nbVertex, 0.0);

            for (int i_Edge = 0; i_Edge < nbEdge; i_Edge++)
            {
                int i_StartVertex = gaussMap.GetEdge(e_VtoG[i_Edge]).AdjacentFace.Index;
                int i_EndVertex = gaussMap.GetEdge(e_VtoG[i_Edge]).PairEdge.AdjacentFace.Index;

                // "X component of Edge" = "X component of EndVertex" - "X component of StartVertex" 
                EfromV[(3 * i_Edge), (3 * i_EndVertex)] = 1.0; EfromV[(3 * i_Edge), (3 * i_StartVertex)] = -1.0;
                // "Y component of Edge" = "Y component of EndVertex" - "Y component of StartVertex" 
                EfromV[(3 * i_Edge) + 1, (3 * i_EndVertex) + 1] = 1.0; EfromV[(3 * i_Edge) + 1, (3 * i_StartVertex) + 1] = -1.0;
                // "Z component of Edge" = "Z component of EndVertex" - "Z component of StartVertex" 
                EfromV[(3 * i_Edge) + 2, (3 * i_EndVertex) + 2] = 1.0; EfromV[(3 * i_Edge) + 2, (3 * i_StartVertex) + 2] = -1.0;
            }

            #endregion

            #region Cross Product with Direction : C

            // With this matrix we want : C * EfromV * V = 0. It expresses to goals :
            // 1. The edge "EfromV * V" of the Voss surface must be parallel to the direction "dir" computed thanks to the Gauss map.
            // 2. The length of the edges "EfromV * V" around a face must be compatible, ie form a closed face. 
            double[,] C = Acc.Matrix.Create(3 * nbEdge, 3 * nbEdge, 0.0);

            for (int i_Edge = 0; i_Edge < nbEdge; i_Edge++)
            {
                Euc.Vector dir = dir_Edges[i_Edge];
                // For the X component of the cross product of dir and edge E[i]
                C[(3 * i_Edge), (3 * i_Edge)] = 0.0; C[(3 * i_Edge), (3 * i_Edge) + 1] = -dir.Z; C[(3 * i_Edge), (3 * i_Edge) + 2] = dir.Y;
                // For the Y component of the cross product of dir and edge E[i]
                C[(3 * i_Edge) + 1, (3 * i_Edge)] = dir.Z; C[(3 * i_Edge) + 1, (3 * i_Edge) + 1] = 0.0; C[(3 * i_Edge) + 1, (3 * i_Edge) + 2] = -dir.X;
                // For the Z component of the cross product of dir and edge E[i]
                C[(3 * i_Edge) + 2, (3 * i_Edge)] = -dir.Y; C[(3 * i_Edge) + 2, (3 * i_Edge) + 1] = dir.X; C[(3 * i_Edge) + 2, (3 * i_Edge) + 2] = 0.0;
            }

            #endregion

            #region Scalar Product with Direction : S

            // With this matrix we want min || S * EfromV * V - L || where L is at target length for each edges.
            // It expresses the fact that the orientation of the Voss edge must be as close as possible from the given orientation of the direction "dir".
            // It avoids having too intricated meshes.

            List<int> g_WWFactors = Enumerable.Repeat(1, gaussMap.EdgeCount / 2).ToList();
            // To allow negative gaussian curvature nets, the warp or weft directions must be inverted.
            if (!Kp)
            {
                int[] g_WarpEdges;
                gaussMap.Quad.WarpAndWeft(out g_WarpEdges, out _);
                for (int i_WarpEdge = 0; i_WarpEdge < g_WarpEdges.Length; i_WarpEdge++)
                {
                    // Problem of conversion between Gauss and voss indices
                    g_WWFactors[g_WarpEdges[i_WarpEdge]] = -1;
                }
            }

            // Computes the matrix
            double[,] S = Acc.Matrix.Create(nbEdge, 3 * nbEdge, 0.0);

            for (int i_Edge = 0; i_Edge < nbEdge; i_Edge++)
            {
                Euc.Vector dir = dir_Edges[i_Edge];
                dir.Unitize();
                // Might invert the direction of the vector when the goal is to compute a surface with negative gaussian curvature
                dir *= g_WWFactors[e_VtoG[i_Edge] / 2];

                // For the scalar product of dir and edge E[i]
                S[(i_Edge), (3 * i_Edge)] = dir.X; S[(i_Edge), (3 * i_Edge) + 1] = dir.Y; S[(i_Edge), (3 * i_Edge) + 2] = dir.Z;
            }

            #endregion


            #region Solver

            // To filter global translations in the modes : We choose to have the vertex "0" at (0,0,0)
            // This results in removing the choice of the first vertex, thus we adapt the matrix EfromV

            double[,] EfromV_f = new double[EfromV.GetLength(0), EfromV.GetLength(1) - 3];
            for (int row = 0; row < EfromV.GetLength(0); row++)
            {
                for (int column = 0; column < EfromV.GetLength(1) - 3; column++)
                {
                    EfromV_f[row, column] = EfromV[row, column + 3];
                }
            }

            // Compute the linear space where V (the coordinate column-vector) can be found.
            // It corresponds to the kernel K of C * EfromV_f).
            double[,] prod = Acc.Matrix.Dot(C, EfromV_f);

            // Get the null space of the matrix "prod"
            var svd = new Acc.Decompositions.SingularValueDecomposition(prod);

            var svd_V = svd.RightSingularVectors;

            var svd_Threshold = Acc.Matrix.GetLength(prod).Max() * Acc.Constants.DoubleEpsilon;
            var E = Acc.Matrix.Find(svd.Diagonal, x => (Double)Math.Abs(x) < svd_Threshold);


            double[,] Kern = Acc.Matrix.GetColumns(svd_V, E);

            /*
                        var svd_Threshold = prod.GetLength().Max() * Acc.Constants.DoubleEpsilon;
                        var E = svd.Diagonal.Find(x => (Double)Math.Abs(x) < svd_Threshold); ;

                        double[,] Kern = svd_V.GetColumns(E);
            */

            // Optimisation to have a "good looking" surface.
            // Here, we want the length of the edges to be close as possible from 1.
            double[] L = Enumerable.Repeat(-1.0, nbEdge).ToArray();


            double[,] SEK = Acc.Matrix.Dot(S, Acc.Matrix.Dot(EfromV_f, Kern));
            double[,] SEK_T = Acc.Matrix.Transpose(SEK);


            double[,] A = Acc.Matrix.Dot(SEK_T, SEK);
            double[] B = Acc.Matrix.Dot(SEK_T, L);


            // Solving : F are the factors for each modes.
            double[] F = Acc.Matrix.Solve(A, B);


            if (F.Length == 0) { return; }
            // If a solution was found, we have to add the vertex "0" at (0,0,0) to the result
            double[] result = Acc.Matrix.Dot(Kern, F);

            // Create the resulting list of vertices
            double[] V = Accord.Math.Vector.Create(result.Length + 3, 0.0);
            V[0] = 0; V[1] = 0; V[2] = 0;
            for (int i = 0; i < result.Length; i++)
            {
                V[i + 3] = result[i];
            }

            #endregion


            #region Generate HeMesh

            for (int i_Vertex = 0; i_Vertex < nbVertex; i_Vertex++)
            {
                vossNet.AddVertex(new Euc.Point(V[(3 * i_Vertex)], V[(3 * i_Vertex) + 1], V[(3 * i_Vertex) + 2]));
            }
            for (int i_Face = 0; i_Face < nbFace; i_Face++)
            {
                List<HeVertex<Euc.Point>> faceVertices = new List<HeVertex<Euc.Point>>();
                foreach (HeEdge<Euc.Point> halfedge in gaussMap.GetVertex(f_VtoG[i_Face]).OutgoingEdges())
                {
                    faceVertices.Add(vossNet.GetVertex(halfedge.AdjacentFace.Index));
                }
                vossNet.AddFace(faceVertices);
            }

            #endregion
        }

        /// <summary> 
        /// Computes Voss net from a Gauss map using a target length for the edges.
        /// </summary>
        /// <param name="gaussMap"> The Gauss map of the Voss net to generate.</param>
        /// <param name="Kp"> Evaluates whether the Gaussian curvature should be positive or not.</param>
        /// <param name="vossNet"> The Voss net generated from the input Gauss map.</param>
        public static void Core_TargetLength_Sparse(HeMesh<Euc.Point> gaussMap, bool Kp, out HeMesh<Euc.Point> vossNet)
        {
            #region Variables Declaration

            vossNet = new HeMesh<Euc.Point>();

            int nbVertex;       // Equal to the number of faces in the Gauss image.
            int nbEdge;         // Equal to the number of internal edges in the Gauss image.
            int nbFace;         // Equal to the number of internal vertices in the Gauss image. 

            List<int> e_VtoG = new List<int>();     // Convert Voss edge indices into Gauss Halfedge indices
            List<int> f_VtoG = new List<int>();     // Convert Voss face indices into Gauss vertex indices

            #endregion

            #region Variables Initialisation

            // Number of Voss vertices
            nbVertex = gaussMap.FaceCount;

            // Number of Voss edges
            for (int i = 0; i < gaussMap.EdgeCount / 2; i++)
            {
                HeEdge<Euc.Point> edge = gaussMap.GetEdge(2 * i);
                if (edge.IsBoundary() || edge.PairEdge.IsBoundary()) { continue; }
                e_VtoG.Add(2 * i);
            }
            nbEdge = e_VtoG.Count;

            //Number of Voss faces
            for (int i = 0; i < gaussMap.VertexCount; i++)
            {
                if (gaussMap.GetVertex(i).IsBoundary()) { continue; }
                f_VtoG.Add(i);
            }
            nbFace = f_VtoG.Count;

            #endregion

            #region Edge Directions 

            Euc.Vector[] dir_Edges = new Euc.Vector[nbEdge];
            for (int i_Edge = 0; i_Edge < nbEdge; i_Edge++)
            {
                Euc.Vector dir = Euc.Vector.CrossProduct((Euc.Vector)gaussMap.GetEdge(e_VtoG[i_Edge]).EndVertex.Position,
                    (Euc.Vector)gaussMap.GetEdge(e_VtoG[i_Edge]).StartVertex.Position);

                dir_Edges[i_Edge] = dir;
            }

            #endregion


            #region Edges from Vertices : EfromV

            // With this matrix we have E = EfromV * V.
            LinAlg.Storage.CoordinateList cl_EfromV_f = new LinAlg.Storage.CoordinateList();

            for (int i_Edge = 0; i_Edge < nbEdge; i_Edge++)
            {
                int i_StartVertex = gaussMap.GetEdge(e_VtoG[i_Edge]).AdjacentFace.Index;
                int i_EndVertex = gaussMap.GetEdge(e_VtoG[i_Edge]).PairEdge.AdjacentFace.Index;

                if (i_StartVertex != 0)
                {
                    i_StartVertex -= 1;
                    // "X component of Edge" = "X component of EndVertex" - "X component of StartVertex" 
                    cl_EfromV_f.Add(-1.0, (3 * i_Edge), (3 * i_StartVertex));
                    // "Y component of Edge" = "Y component of EndVertex" - "Y component of StartVertex" 
                    cl_EfromV_f.Add(-1.0, (3 * i_Edge) + 1, (3 * i_StartVertex) + 1);
                    // "Z component of Edge" = "Z component of EndVertex" - "Z component of StartVertex" 
                    cl_EfromV_f.Add(-1.0, (3 * i_Edge) + 2, (3 * i_StartVertex) + 2);
                }
                if (i_EndVertex != 0)
                {
                    i_EndVertex -= 1;
                    // "X component of Edge" = "X component of EndVertex" - "X component of StartVertex" 
                    cl_EfromV_f.Add(1.0, (3 * i_Edge), (3 * i_EndVertex));
                    // "Y component of Edge" = "Y component of EndVertex" - "Y component of StartVertex" 
                    cl_EfromV_f.Add(1.0, (3 * i_Edge) + 1, (3 * i_EndVertex) + 1);
                    // "Z component of Edge" = "Z component of EndVertex" - "Z component of StartVertex" 
                    cl_EfromV_f.Add(1.0, (3 * i_Edge) + 2, (3 * i_EndVertex) + 2);
                }
            }

            LinAlg.SparseMatrix EfromV_f = cl_EfromV_f.ToSparseMatrix(3 * nbEdge, 3 * (nbVertex - 1));

            #endregion

            #region Cross Product with Direction : C

            // With this matrix we want : C * EfromV * V = 0. It expresses to goals :
            // 1. The edge "EfromV * V" of the Voss surface must be parallel to the direction "dir" computed thanks to the Gauss map.
            // 2. The length of the edges "EfromV * V" around a face must be compatible, ie form a closed face. 
            LinAlg.Storage.CoordinateList cl_C = new LinAlg.Storage.CoordinateList();

            for (int i_Edge = 0; i_Edge < nbEdge; i_Edge++)
            {
                Euc.Vector dir = dir_Edges[i_Edge];
                // For the X component of the cross product of dir and edge E[i]
                cl_C.Add(-dir.Z, (3 * i_Edge), (3 * i_Edge) + 1); cl_C.Add(dir.Y, (3 * i_Edge), (3 * i_Edge) + 2);
                // For the Y component of the cross product of dir and edge E[i]
                cl_C.Add(dir.Z, (3 * i_Edge) + 1, (3 * i_Edge)); cl_C.Add(-dir.X, (3 * i_Edge) + 1, (3 * i_Edge) + 2);
                // For the Z component of the cross product of dir and edge E[i]
                cl_C.Add(-dir.Y, (3 * i_Edge) + 2, (3 * i_Edge)); cl_C.Add(dir.X, (3 * i_Edge) + 2, (3 * i_Edge) + 1);
            }


            LinAlg.SparseMatrix C = cl_C.ToSparseMatrix(3 * nbEdge, 3 * nbEdge);

            #endregion

            #region Scalar Product with Direction : S

            // With this matrix we want min || S * EfromV * V - L || where L is at target length for each edges.
            // It expresses the fact that the orientation of the Voss edge must be as close as possible from the given orientation of the direction "dir".
            // It avoids having too intricated meshes.

            List<int> g_WWFactors = Enumerable.Repeat(1, gaussMap.EdgeCount / 2).ToList();
            // To allow negative gaussian curvature nets, the warp or weft directions must be inverted.
            if (!Kp)
            {
                int[] g_WarpEdges;
                gaussMap.Quad.WarpAndWeft(out g_WarpEdges, out _);
                for (int i_WarpEdge = 0; i_WarpEdge < g_WarpEdges.Length; i_WarpEdge++)
                {
                    // Problem of conversion between Gauss and voss indices
                    g_WWFactors[g_WarpEdges[i_WarpEdge]] = -1;
                }
            }

            // Computes the matrix
            LinAlg.Storage.CoordinateList cl_S = new LinAlg.Storage.CoordinateList();

            for (int i_Edge = 0; i_Edge < nbEdge; i_Edge++)
            {
                Euc.Vector dir = dir_Edges[i_Edge];
                dir.Unitize();
                // Might invert the direction of the vector when the goal is to compute a surface with negative gaussian curvature
                dir *= g_WWFactors[e_VtoG[i_Edge] / 2];

                // For the scalar product of dir and edge E[i]
                cl_S.Add(dir.X, (i_Edge), (3 * i_Edge)); cl_S.Add(dir.Y, (i_Edge), (3 * i_Edge) + 1); cl_S.Add(dir.Z, (i_Edge), (3 * i_Edge) + 2);
            }


            LinAlg.SparseMatrix S = cl_S.ToSparseMatrix(nbEdge, 3 * nbEdge);

            #endregion


            #region Solver

            // Compute the linear space where V (the coordinate column-vector) can be found.
            // It corresponds to the kernel K of C * EfromV_f).

            LinAlg.SparseMatrix prod = LinAlg.SparseMatrix.Multiply(C, EfromV_f);
            LinAlg.Vector[] kern = LinAlg.SparseMatrix.Kernel(prod);

            LinAlg.Storage.CoordinateList cl_K = new LinAlg.Storage.CoordinateList();
            for (int i_Vect = 0; i_Vect < kern.Length; i_Vect++)
            {
                int dim = kern[i_Vect].Size;
                for (int i_Coord = 0; i_Coord < dim; i_Coord++)
                {
                    if (Math.Abs(kern[i_Vect][i_Coord]) > 1e-15)                 // BEWARE OF THIS CONDITION
                    {
                        cl_K.Add(kern[i_Vect][i_Coord], i_Coord, i_Vect);
                    }
                }
            }

            LinAlg.SparseMatrix K = cl_K.ToSparseMatrix(EfromV_f.ColumnCount, kern.Length);

            // Optimisation to have a "good looking" surface.
            // Here, we want the length of the edges to be close as possible from 1.

            LinAlg.Vector L = new LinAlg.Vector(nbEdge, -1.0);

            LinAlg.SparseMatrix Temp = LinAlg.SparseMatrix.Multiply(EfromV_f, K);
            LinAlg.SparseMatrix SEK = LinAlg.SparseMatrix.Multiply(S, Temp);

            LinAlg.SparseMatrix A = LinAlg.SparseMatrix.Square(SEK);
            LinAlg.Vector B = LinAlg.SparseMatrix.TransposeMultiply(SEK, L);

            // Solving : F are the factors for each modes.
            LinAlg.Vector F = A.SolveQR(B);

            if (F.Size == 0) { return; }
            // If a solution was found, we have to add the vertex "0" at (0,0,0) to the result
            LinAlg.Vector V = LinAlg.SparseMatrix.Multiply(K, F);

            #endregion


            #region Generate HeMesh

            for (int i_Vertex = 0; i_Vertex < nbVertex; i_Vertex++)
            {
                vossNet.AddVertex(new Euc.Point(V[(3 * i_Vertex)], V[(3 * i_Vertex) + 1], V[(3 * i_Vertex) + 2]));
            }
            for (int i_Face = 0; i_Face < nbFace; i_Face++)
            {
                List<HeVertex<Euc.Point>> faceVertices = new List<HeVertex<Euc.Point>>();
                foreach (HeEdge<Euc.Point> halfedge in gaussMap.GetVertex(f_VtoG[i_Face]).OutgoingEdges())
                {
                    faceVertices.Add(vossNet.GetVertex(halfedge.AdjacentFace.Index));
                }
                vossNet.AddFace(faceVertices);
            }

            #endregion
        }

        /******************** DEveloppable Voss Net ********************/

        /// <summary>
        /// Computes a Developable/Voss net form its Gauss image on the unity sphere.
        /// </summary>
        /// <param name="gaussMap"> Gauss image of the Developpable net.</param>
        /// <param name="faceCount"> Number of faces in the other direction.</param>
        /// <param name="devNet"> Developpable net.</param>
        [Obsolete]
        public static void GaussMapToDevelopableNet_Accord(List<Euc.Point> gaussMap, int faceCount, out HeMesh<Euc.Point> devNet)
        {
            if (DateTime.Now > new DateTime(2020, 08, 31)) { throw new UnauthorizedAccessException("The license of this component is revoqued."); }

            #region On DevNet

            devNet = new HeMesh<Euc.Point>();
            // Face
            int nb_U_Face = faceCount;
            int nb_V_Face = gaussMap.Count;
            // Edge
            int nbEdge_U = nb_U_Face * (nb_V_Face + 1);
            int nbEdge = (nb_U_Face * (nb_V_Face + 1)) + (nb_V_Face * (nb_U_Face + 1));
            // Vertex
            int nb_Vertex = (nb_U_Face + 1) * (nb_V_Face + 1);

            Euc.Vector[] u_Directions = new Euc.Vector[nb_V_Face + 1];
            Euc.Vector[] v_Directions = new Euc.Vector[nb_V_Face];

            #endregion


            #region Compute Edge Directions

            // First "U" direction
            Euc.Vector dir_First = Euc.Vector.CrossProduct((Euc.Vector)gaussMap[0], (Euc.Vector)gaussMap[1]);
            dir_First.Unitize();
            u_Directions[0] = new Euc.Vector(dir_First);
            // "U" directions
            for (int i = 0; i < nb_V_Face - 1; i++)
            {
                Euc.Vector dir = Euc.Vector.CrossProduct((Euc.Vector)gaussMap[i], (Euc.Vector)gaussMap[i + 1]);
                dir.Unitize();
                u_Directions[i + 1] = new Euc.Vector(dir);
            }
            // Last "U" direction
            u_Directions[nb_V_Face] = new Euc.Vector(u_Directions[nb_V_Face - 1]);

            // "V" directions
            for (int i = 0; i < nb_V_Face; i++)
            {
                v_Directions[i] = Euc.Vector.CrossProduct((Euc.Vector) gaussMap[i], u_Directions[i]);
            }

            #endregion

            #region Edges from Vertices : EfromV

            // With this matrix we have E = EfromV * V.
            double[,] EfromV = Acc.Matrix.Create(3 * nbEdge, 3 * nb_Vertex, 0.0);

            // For edges in the "U" directions
            for (int row = 0; row < nb_V_Face + 1; row++)
            {
                for (int rank = 0; rank < nb_U_Face; rank++)
                {
                    // Get start and end vertex indices
                    int i_StartVertex = row * (nb_U_Face + 1) + rank;
                    int i_EndVertex = i_StartVertex + 1;

                    int i_Edge = (row * nb_U_Face) + rank;
                    // "X component of Edge" = "X component of EndVertex" - "X component of StartVertex" 
                    EfromV[(3 * i_Edge), (3 * i_EndVertex)] = 1.0; EfromV[(3 * i_Edge), (3 * i_StartVertex)] = -1.0;
                    // "Y component of Edge" = "Y component of EndVertex" - "Y component of StartVertex" 
                    EfromV[(3 * i_Edge) + 1, (3 * i_EndVertex) + 1] = 1.0; EfromV[(3 * i_Edge) + 1, (3 * i_StartVertex) + 1] = -1.0;
                    // "Z component of Edge" = "Z component of EndVertex" - "Z component of StartVertex" 
                    EfromV[(3 * i_Edge) + 2, (3 * i_EndVertex) + 2] = 1.0; EfromV[(3 * i_Edge) + 2, (3 * i_StartVertex) + 2] = -1.0;
                }
            }
            // For edges in the "V" directions
            for (int row = 0; row < nb_V_Face; row++)
            {
                for (int rank = 0; rank < nb_U_Face + 1; rank++)
                {
                    // Get start and end vertex indices
                    int i_StartVertex = (row * (nb_U_Face + 1)) + rank;
                    int i_EndVertex = i_StartVertex + (nb_U_Face + 1);

                    int i_Edge = nbEdge_U + i_StartVertex;
                    // "X component of Edge" = "X component of EndVertex" - "X component of StartVertex" 
                    EfromV[(3 * i_Edge), (3 * i_EndVertex)] = 1.0; EfromV[(3 * i_Edge), (3 * i_StartVertex)] = -1.0;
                    // "Y component of Edge" = "Y component of EndVertex" - "Y component of StartVertex" 
                    EfromV[(3 * i_Edge) + 1, (3 * i_EndVertex) + 1] = 1.0; EfromV[(3 * i_Edge) + 1, (3 * i_StartVertex) + 1] = -1.0;
                    // "Z component of Edge" = "Z component of EndVertex" - "Z component of StartVertex" 
                    EfromV[(3 * i_Edge) + 2, (3 * i_EndVertex) + 2] = 1.0; EfromV[(3 * i_Edge) + 2, (3 * i_StartVertex) + 2] = -1.0;
                }
            }

            #endregion

            #region Cross Product with Direction : C

            // With this matrix we want : C * EfromV * V = 0. It expresses to goals :
            // 1. The edge "EfromV * V" of the Voss surface must be parallel to the direction "dir" computed thanks to the Gauss map.
            // 2. The length of the edges "EfromV * V" around a face must be compatible, ie form a closed face. 
            double[,] C = Acc.Matrix.Create(3 * nbEdge, 3 * nbEdge, 0.0);

            // For edges in the "U" directions
            for (int row = 0; row < nb_V_Face + 1; row++)
            {
                for (int rank = 0; rank < nb_U_Face; rank++)
                {
                    // Get the direction for the edge
                    Euc.Vector dir = u_Directions[row];

                    int i_Edge = (row * nb_U_Face) + rank;
                    // For the X component of the cross product of dir and edge E[i]
                    C[(3 * i_Edge), (3 * i_Edge)] = 0; C[(3 * i_Edge), (3 * i_Edge) + 1] = -dir.Z; C[(3 * i_Edge), (3 * i_Edge) + 2] = dir.Y;
                    // For the Y component of the cross product of dir and edge E[i]
                    C[(3 * i_Edge) + 1, (3 * i_Edge)] = dir.Z; C[(3 * i_Edge) + 1, (3 * i_Edge) + 1] = 0; C[(3 * i_Edge) + 1, (3 * i_Edge) + 2] = -dir.X;
                    // For the Z component of the cross product of dir and edge E[i]
                    C[(3 * i_Edge) + 2, (3 * i_Edge)] = -dir.Y; C[(3 * i_Edge) + 2, (3 * i_Edge) + 1] = dir.X; C[(3 * i_Edge) + 2, (3 * i_Edge) + 2] = 0;
                }
            }
            // For edges in the "V" directions
            for (int row = 0; row < nb_V_Face; row++)
            {
                for (int rank = 0; rank < nb_U_Face + 1; rank++)
                {
                    // Get the direction for the edge
                    Euc.Vector dir = v_Directions[row];

                    int i_Edge = nbEdge_U + (row * (nb_U_Face + 1)) + rank;
                    // For the X component of the cross product of dir and edge E[i]
                    C[(3 * i_Edge), (3 * i_Edge)] = 0; C[(3 * i_Edge), (3 * i_Edge) + 1] = -dir.Z; C[(3 * i_Edge), (3 * i_Edge) + 2] = dir.Y;
                    // For the Y component of the cross product of dir and edge E[i]
                    C[(3 * i_Edge) + 1, (3 * i_Edge)] = dir.Z; C[(3 * i_Edge) + 1, (3 * i_Edge) + 1] = 0; C[(3 * i_Edge) + 1, (3 * i_Edge) + 2] = -dir.X;
                    // For the Z component of the cross product of dir and edge E[i]
                    C[(3 * i_Edge) + 2, (3 * i_Edge)] = -dir.Y; C[(3 * i_Edge) + 2, (3 * i_Edge) + 1] = dir.X; C[(3 * i_Edge) + 2, (3 * i_Edge) + 2] = 0;
                }
            }

            #endregion

            #region Scalar Product with Direction : S

            // With this matrix we want min || S * EfromV * V - L || where L is at target length for each edges.
            // It expresses the fact that the orientation of the Voss edge must be as close as possible from the given orientation of the direction "dir".
            // It avoids having too intricated meshes.

            // Computes the matrix
            double[,] S = Acc.Matrix.Create(nbEdge, 3 * nbEdge, 0.0);

            // For edges in the "U" directions
            for (int row = 0; row < nb_V_Face + 1; row++)
            {
                for (int rank = 0; rank < nb_U_Face; rank++)
                {
                    // Get the direction for the edge
                    Euc.Vector dir = u_Directions[row];

                    int i_Edge = (row * nb_U_Face) + rank;
                    // For the scalar product of dir and edge E[i]
                    S[(i_Edge), (3 * i_Edge)] = dir.X; S[(i_Edge), (3 * i_Edge) + 1] = dir.Y; S[(i_Edge), (3 * i_Edge) + 2] = dir.Z;
                }
            }
            // For edges in the "V" directions
            for (int row = 0; row < nb_V_Face; row++)
            {
                for (int rank = 0; rank < nb_U_Face + 1; rank++)
                {
                    // Get the direction for the edge
                    Euc.Vector dir = v_Directions[row];

                    int i_Edge = nbEdge_U + (row * (nb_U_Face + 1)) + rank;
                    // For the scalar product of dir and edge E[i]
                    S[(i_Edge), (3 * i_Edge)] = dir.X; S[(i_Edge), (3 * i_Edge) + 1] = dir.Y; S[(i_Edge), (3 * i_Edge) + 2] = dir.Z;
                }
            }

            #endregion


            #region Solver

            // To filter global translations in the modes : We choose to have the vertex "0" at (0,0,0)
            // This results in removing the choice of the first vertex, thus we adapt the matrix EfromV
            double[,] EfromV_f = Acc.Matrix.Create(EfromV.GetLength(0), EfromV.GetLength(1) - 3, 0.0);
            for (int row = 0; row < EfromV.GetLength(0); row++)
            {
                for (int column = 0; column < EfromV.GetLength(1) - 3; column++)
                {
                    EfromV_f[row, column] = EfromV[row, column + 3];
                }
            }


            // Compute the linear space where V (the oordinate column-vector) can be found.
            // It corresponds to the kernel K of C * EfromV_f).
            double[,] C_E_f = Acc.Matrix.Dot(C, EfromV_f);

            // Get the null space of the matrix "prod"
            var svd = new Accord.Math.Decompositions.SingularValueDecomposition(C_E_f);

            var svd_V = svd.RightSingularVectors;

            var svd_Threshold = Acc.Matrix.GetLength(C_E_f).Max() * Acc.Constants.DoubleEpsilon;
            var E = Acc.Matrix.Find(svd.Diagonal, x => (Double)Math.Abs(x) < svd_Threshold); ;

            double[,] Kern = Acc.Matrix.GetColumns(svd_V, E);

            // Optimisation to have a "good looking" surface.
            // Here, we want the length of the edges to be close as possible from 1.
            double[] L = Enumerable.Repeat(-1.0, nbEdge).ToArray();

            double[,] S_E_K = Acc.Matrix.Dot(S, Acc.Matrix.Dot(EfromV_f, Kern));


            double[,] A = Acc.Matrix.Dot(Acc.Matrix.Transpose(S_E_K), S_E_K);
            double[] B = Acc.Matrix.Dot(Acc.Matrix.Transpose(S_E_K), L);


            // Solving : F are the factors for each modes.
            double[] F = Acc.Matrix.Solve(A,B);


            if (F.Length == 0) { return; }
            // If a solution was found, we have to add the vertex "0" at (0,0,0) to the result
            double[] result = Acc.Matrix.Dot(Kern, F);

            // Create the resulting list of vertices
            double[] V = new double[result.Length + 3];
            V[0] = 0; V[1] = 0; V[2] = 0;
            for (int i = 0; i < result.Length; i++)
            {
                V[i + 3] = result[i];
            }

            #endregion


            #region Generate HeMesh

            for (int i_Vertex = 0; i_Vertex < nb_Vertex; i_Vertex++)
            {
                devNet.AddVertex(new Euc.Point(V[(3 * i_Vertex)], V[(3 * i_Vertex) + 1], V[(3 * i_Vertex) + 2]));
            }
            for (int i_Face = 0; i_Face < nb_U_Face * nb_V_Face; i_Face++)
            {
                int row = i_Face / nb_U_Face;
                int rank = i_Face % nb_U_Face;

                int v1 = row * (nb_U_Face + 1) + rank;
                int v2 = v1 + 1;
                int v3 = (row + 1) * (nb_U_Face + 1) + rank + 1;
                int v4 = v3 - 1;


                devNet.AddFace(devNet.GetVertex(v1), devNet.GetVertex(v2), devNet.GetVertex(v3), devNet.GetVertex(v4));
            }

            #endregion

        }
    }
}
