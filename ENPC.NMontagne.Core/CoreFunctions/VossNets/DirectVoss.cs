using System;
using System.Collections.Generic;

using Euc = ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;
using ENPC.Numerics;


namespace ENPC.NMontagne.Core.CoreFunctions.VossNets
{
    /// <summary>
    /// Class containing methods to generate Voss nets from boundary curves.
    /// </summary>
    public static class DirectVoss
    {

        /// <summary>
        /// Computes a Voss net from two primal conditions (two geodesic curves).
        /// </summary>
        /// <param name="NodeNormals_U"> The projection on the unit sphere of the first curve's frenet frame normals.</param>
        /// <param name="NodeNormals_V"> The projection on the unit sphere of the second curve's frenet frame normals.</param>
        /// <param name="Edges_U"> The segments of the first curve.</param>
        /// <param name="Edges_V"> The segments of the second curve.</param>
        /// <param name="gaussMap"> The Gauss map obtained from the primal conditions. </param>
        /// <param name="vossNet"> The Voss net obtained from the primal conditions.</param>
        [Obsolete]
        public static void Core_FromTwoPrimal(List<Euc.Point> NodeNormals_U, List<Euc.Point> NodeNormals_V, List<Euc.Vector> Edges_U, List<Euc.Vector> Edges_V,
            out HeMesh<Euc.Point> gaussMap, out HeMesh<Euc.Point> vossNet)
        {
            #region Verifications

            /********** Verifications **********/
            if (NodeNormals_U[0].DistanceTo(NodeNormals_V[0]) > Settings._absolutePrecision)              // NodeNormal U-V coincidence
            {
                throw new ArgumentException("The Node Normals U and V must coincide at start.");
            }
            if (Edges_U.Count + 1 != NodeNormals_U.Count || Edges_V.Count + 1 != NodeNormals_V.Count)      // Number of U-V elements 
            {
                throw new ArgumentException("The number of edges must be one less than the number of normals in each directions.");
            }

            /********** Verifes NodeNormal U-V Orientation **********/

            Euc.Vector temp_U = Euc.Vector.CrossProduct((Euc.Vector) NodeNormals_U[0], (Euc.Vector) NodeNormals_U[1]);
            Euc.Vector temp_V = Euc.Vector.CrossProduct((Euc.Vector) NodeNormals_V[0], (Euc.Vector) NodeNormals_V[1]);

            Euc.Vector temp_N = Euc.Vector.CrossProduct(temp_U, temp_V);

            if (Euc.Vector.DotProduct(temp_N, (Euc.Vector) NodeNormals_U[0]) < 0)
            {
                List<Euc.Point> tmp_PointList = new List<Euc.Point>(NodeNormals_U);
                NodeNormals_U = NodeNormals_V;
                NodeNormals_V = tmp_PointList;

                List<Euc.Vector> tmp_VectorList = new List<Euc.Vector>(Edges_U);
                Edges_U = Edges_V;
                Edges_V = tmp_VectorList;
            }

            /********** Arrange Edges Orientation **********/

            double[] EdgesOri_U = new double[Edges_U.Count];
            for (int i = 0; i < Edges_U.Count; i++)
            {
                Euc.Vector expectedOri = (Euc.Vector) (NodeNormals_U[i + 1] - NodeNormals_U[i]);
                if (Euc.Vector.DotProduct(expectedOri, Edges_U[i]) < 0) { EdgesOri_U[i] = -1 / Edges_U[i].Length(); }
                else { EdgesOri_U[i] = 1 / Edges_U[i].Length(); }
            }

            double[] EdgesOri_V = new double[Edges_V.Count];
            for (int i = 0; i < Edges_V.Count; i++)
            {
                Euc.Vector expectedOri = (Euc.Vector) (NodeNormals_V[i + 1] - NodeNormals_V[i]);
                if (Euc.Vector.DotProduct(expectedOri, Edges_V[i]) < 0) { EdgesOri_V[i] = -1 / Edges_V[i].Length(); }
                else { EdgesOri_V[i] = 1 / Edges_V[i].Length(); }
            }

            #endregion

            #region Initialisation

            /********** Computation of the first FaceNormal **********/

            Euc.Vector FaceNormal_1_1 = Euc.Vector.CrossProduct(EdgesOri_U[0] * Edges_U[0], EdgesOri_V[0] * Edges_V[0]);
            FaceNormal_1_1.Unitize();

            // Get the first face normal
            ((Euc.Vector) NodeNormals_U[0]).Unitize();
            Quaternion qAxis = Quaternion.Unit(((Euc.Vector)NodeNormals_U[0]), Math.PI);
            Quaternion qAxis_1 = qAxis.Inverse();
            Quaternion qVector = new Quaternion(0, FaceNormal_1_1.X, FaceNormal_1_1.Y, FaceNormal_1_1.Z);
            qVector = qAxis * (qVector * qAxis_1);

            Euc.Vector FaceNormal_0_0 = new Euc.Vector(qVector.I, qVector.J, qVector.K);

            List<Euc.Vector> FaceNormals_U = new List<Euc.Vector>(); FaceNormals_U.Add(FaceNormal_0_0);
            List<Euc.Vector> FaceNormals_V = new List<Euc.Vector>(); FaceNormals_V.Add(FaceNormal_0_0);

            /********** Computation the first dihedral angles **********/

            Euc.Vector tangent = Euc.Vector.CrossProduct(FaceNormal_1_1, Euc.Vector.CrossProduct(FaceNormal_0_0, FaceNormal_1_1));
            tangent.Unitize();
            double a = Math.Acos(Euc.Vector.DotProduct(FaceNormal_0_0, FaceNormal_1_1));

            Euc.Vector forBeta = Euc.Vector.CrossProduct(EdgesOri_U[0] * Edges_U[0], FaceNormal_1_1);
            forBeta.Unitize();
            double beta = Math.Acos(Euc.Vector.DotProduct(forBeta, tangent));

            Euc.Vector forGamma = Euc.Vector.CrossProduct(FaceNormal_1_1, EdgesOri_V[0] * Edges_V[0]);
            forGamma.Unitize();
            double gamma = Math.Acos(Euc.Vector.DotProduct(forGamma, tangent));

            double alpha = Math.Acos(-(Math.Cos(beta) * Math.Cos(gamma)) + (Math.Sin(beta) * Math.Sin(gamma) * Math.Cos(a)));

            double mu = Math.Asin((Math.Sin(a) * Math.Sin(gamma)) / Math.Sin(alpha));
            double nu = Math.Asin((Math.Sin(a) * Math.Sin(beta)) / Math.Sin(alpha));

            #endregion

            #region Compute Chebyshev Boundaries (Iteration)

            /********** For the FaceNormals_U **********/

            for (int i = 0; i < NodeNormals_U.Count - 1; i++)
            {
                qAxis = Quaternion.Unit(((Euc.Vector) NodeNormals_U[i]), Math.PI);
                qAxis_1 = qAxis.Inverse();
                qVector = new Quaternion(0, FaceNormals_U[i].X, FaceNormals_U[i].Y, FaceNormals_U[i].Z);
                qVector = qAxis * qVector * qAxis_1;

                qAxis = Quaternion.Unit(EdgesOri_U[i] * Edges_U[i], mu);
                qAxis_1 = qAxis.Inverse();
                qVector = qAxis * qVector * qAxis_1;

                FaceNormals_U.Add(new Euc.Vector(qVector.I, qVector.J, qVector.K));
            }

            // For the last FaceNormals_U
            qVector = qAxis_1 * qVector * qAxis;

            qAxis = Quaternion.Unit(((Euc.Vector) NodeNormals_U[NodeNormals_U.Count - 1]), Math.PI);
            qAxis_1 = qAxis.Inverse();
            qVector = qAxis * qVector * qAxis_1;

            FaceNormals_U.Add(new Euc.Vector(qVector.I, qVector.J, qVector.K));


            /********** For the FaceNormals_V **********/

            for (int i = 0; i < NodeNormals_V.Count - 1; i++)
            {
                qAxis = Quaternion.Unit(((Euc.Vector) NodeNormals_V[i]), Math.PI);
                qAxis_1 = qAxis.Inverse();
                qVector = new Quaternion(0, FaceNormals_V[i].X, FaceNormals_V[i].Y, FaceNormals_V[i].Z);
                qVector = qAxis * qVector * qAxis_1;

                qAxis = Quaternion.Unit(EdgesOri_V[i] * Edges_V[i], -nu);
                qAxis_1 = qAxis.Inverse();
                qVector = qAxis * qVector * qAxis_1;

                FaceNormals_V.Add(new Euc.Vector(qVector.I, qVector.J, qVector.K));
            }

            // For the last FaceNormals_V
            qVector = qAxis_1 * qVector * qAxis;

            qAxis = Quaternion.Unit(((Euc.Vector) NodeNormals_V[NodeNormals_V.Count - 1]), Math.PI);
            qAxis_1 = qAxis.Inverse();
            qVector = qAxis * qVector * qAxis_1;

            FaceNormals_V.Add(new Euc.Vector(qVector.I, qVector.J, qVector.K));

            #endregion

            #region Compute Chebyhev From Primals

            Euc.Vector[,] FaceNormals = new Euc.Vector[FaceNormals_U.Count, FaceNormals_V.Count];

            for (int i = 0; i < FaceNormals_U.Count; i++) { FaceNormals[i, 0] = FaceNormals_U[i]; }
            for (int j = 0; j < FaceNormals_V.Count; j++) { FaceNormals[0, j] = FaceNormals_V[j]; }

            int maxSum = FaceNormals_V.Count + FaceNormals_U.Count;
            for (int sum = 2; sum < maxSum; sum++)
            {
                for (int i = 1; i < sum; i++)
                {
                    int j = sum - i;
                    if (!(i < FaceNormals_U.Count) || !(j < FaceNormals_V.Count)) { continue; }
                    // Rotation Quaternions
                    Euc.Vector axis = (FaceNormals[i - 1, j] + FaceNormals[i, j - 1]);
                    axis.Unitize();
                    double angle = Math.PI;
                    qAxis = new Quaternion(Math.Cos(angle / 2),
                      Math.Sin(angle / 2) * axis.X, Math.Sin(angle / 2) * axis.Y, Math.Sin(angle / 2) * axis.Z);
                    qAxis_1 = qAxis.Inverse();
                    // Vector Quaternion
                    qVector = new Quaternion(0, FaceNormals[i - 1, j - 1].X,
                      FaceNormals[i - 1, j - 1].Y, FaceNormals[i - 1, j - 1].Z);
                    // Compute rotation about "axis" of "angle" of "vector"
                    Quaternion qResult = qAxis * qVector * qAxis_1;
                    FaceNormals[i, j] = new Euc.Vector(qResult.I, qResult.J, qResult.K);
                }
            }

            #endregion

            #region Voss from Chebyshev and Geodesics ("Linear" but Iterative)

            Euc.Point[,] vossVertices = new Euc.Point[Edges_U.Count + 1, Edges_V.Count + 1];

            // Initialisation
            vossVertices[0, 0] = new Euc.Point(0, 0, 0);
            for (int i = 0; i < Edges_U.Count; i++)
            {
                vossVertices[i + 1, 0] = vossVertices[i, 0] + Edges_U[i];
            }
            for (int j = 0; j < Edges_V.Count; j++)
            {
                vossVertices[0, j + 1] = vossVertices[0, j] + Edges_V[j];
            }

            // Fill array of vossVertices
            maxSum = vossVertices.GetLength(0) + vossVertices.GetLength(1);
            for (int sum = 2; sum < maxSum; sum++)
            {
                for (int i = 1; i < sum; i++)
                {
                    int j = sum - i;
                    // Verifications
                    if (!(i < Edges_U.Count + 1) || !(j < Edges_V.Count + 1)) { continue; }

                    // Computes Edges
                    Euc.Vector dirB, dirC;

                    dirB = Euc.Vector.CrossProduct(FaceNormals[i + 1, j], FaceNormals[i, j]);
                    dirB.Unitize();
                    //if (Vector3d.DotProduct(dirB, Edges_V[j - 1]) < 0) { dirB = -dirB; }


                    dirC = Euc.Vector.CrossProduct(FaceNormals[i, j + 1], FaceNormals[i, j]);
                    dirC.Unitize();
                    //if (Vector3d.DotProduct(dirC, Edges_U[i - 1]) > 0) { dirC = - dirC; }

                    Euc.Vector diagA = (Euc.Vector)(vossVertices[i, j - 1] - vossVertices[i - 1, j]);
                    double a2 = diagA.Length();
                    diagA.Unitize();

                    beta = Math.PI - Math.Acos(Euc.Vector.DotProduct(diagA, dirC));
                    gamma = Math.PI - Math.Acos(Euc.Vector.DotProduct(diagA, dirB));

                    alpha = Math.PI - Math.Acos(Euc.Vector.DotProduct(dirB, dirC));

                    double b = (a2 * Math.Sin(beta)) / Math.Sin(alpha);
                    double c = (a2 * Math.Sin(gamma)) / Math.Sin(alpha);

                    Euc.Point resultB = vossVertices[i, j - 1] + (b * dirB);
                    Euc.Point resultC = vossVertices[i - 1, j] + (c * dirC);

                    if (resultB.DistanceTo(resultC) < 1) { /*Do nothing*/ }
                    else
                    {
                        resultB = vossVertices[i, j - 1] - (b * dirB);
                        resultC = vossVertices[i - 1, j] + (c * dirC);
                        if (resultB.DistanceTo(resultC) < 1) { /*Do nothing*/ }
                        else
                        {
                            resultB = vossVertices[i, j - 1] + (b * dirB);
                            resultC = vossVertices[i - 1, j] - (c * dirC);
                            double dist = resultB.DistanceTo(resultC);
                            if (resultB.DistanceTo(resultC) < 1) { /*Do nothing*/ }
                            else
                            {
                                resultB = vossVertices[i, j - 1] - (b * dirB);
                                resultC = vossVertices[i - 1, j] - (c * dirC);
                                if (resultB.DistanceTo(resultC) < 1) { /*Do nothing*/ }
                            }
                        }
                    }

                    vossVertices[i, j] = resultB;
                }
            }

            #endregion

            #region Generate Meshes

            gaussMap = new HeMesh<Euc.Point>();

            int nb_U = Edges_U.Count + 2;
            int nb_V = Edges_V.Count + 2;

            for (int i = 0; i < nb_U; i++)
            {
                for (int j = 0; j < nb_V; j++)
                {
                    gaussMap.AddVertex((Euc.Point)FaceNormals[i, j]);
                }
            }

            for (int i = 0; i < nb_U - 1; i++)
            {
                for (int j = 0; j < nb_V -1; j++)
                {
                    gaussMap.AddFace(gaussMap.GetVertex(i + (j * nb_U)), gaussMap.GetVertex((i+1) + (j * nb_U)), gaussMap.GetVertex((i+1) + ((j+1) * nb_U)), gaussMap.GetVertex(i + ((j+1) * nb_U)));
                }
            }



            nb_U = Edges_U.Count + 1;
            nb_V = Edges_V.Count + 1;

            vossNet = new HeMesh<Euc.Point>();
            for (int i = 0; i < nb_U; i++)
            {
                for (int j = 0; j < nb_V; j++)
                {
                    vossNet.AddVertex(vossVertices[i, j]);
                }
            }

            for (int i = 0; i < nb_U - 1; i++)
            {
                for (int j = 0; j < nb_V - 1; j++)
                {
                    gaussMap.AddFace(gaussMap.GetVertex(i + (j * nb_U)), gaussMap.GetVertex((i + 1) + (j * nb_U)), gaussMap.GetVertex((i + 1) + ((j + 1) * nb_U)), gaussMap.GetVertex(i + ((j + 1) * nb_U)));
                }
            }

            #endregion
        }
    }
}
