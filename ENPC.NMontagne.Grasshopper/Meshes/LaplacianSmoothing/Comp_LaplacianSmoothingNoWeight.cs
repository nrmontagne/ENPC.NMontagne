using System;
using System.Collections.Generic;

using Euc = ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;

using GH_K = Grasshopper.Kernel;

using Param_Euc = ENPC.McNeel.Grasshopper.Parameters.Euclidean;

using ENPC.NMontagne.Core.CoreFunctions.Meshes;
using Rhino.Geometry;


namespace ENPC.NMontagne.Grasshopper.Meshes
{
    /// <summary>
    /// A grasshopper component computing the laplacian smoothing of a net.
    /// </summary>
    public class Comp_LaplacianSmoothingNoWeight : GH_K.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initializes a new instance of the <see cref="Comp_LaplacianSmoothingNoWeight"/> class.
        /// </summary>
        public Comp_LaplacianSmoothingNoWeight()
          : base("Laplacian Smoothing", "Lap.",
              "Computes the Laplacian smothing algorithm on a net.",
              LibrarySettings.Otter.Name, LibrarySettings.Otter.Meshes.Name)
        {
        }

        #endregion

        #region Override : GH_Component

        /// <inheritdoc cref="GH_K.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_K.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Param_Euc.Param_HeMesh(), "Net", "M", "The net to operate on.", GH_K.GH_ParamAccess.item);
            pManager.AddIntegerParameter("Iteration", "I","Number of iteration for the algorithm", GH_K.GH_ParamAccess.list);
            pManager.AddIntegerParameter("Boundary Conditions", "C", "The boundary condition : /n 0: free edges; 1: fixed boundary;", GH_K.GH_ParamAccess.list);
        }

        /// <inheritdoc cref="GH_K.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_K.GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Param_Euc.Param_HeMesh(), "Smoothed Net", "N", "The net resulting from the laplacian smoothing.", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.SolveInstance(GH_K.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_K.IGH_DataAccess DA)
        {
            // Declaration, initialization, and instanciation of input variables
            HeMesh<Euc.Point> mesh = new HeMesh<Euc.Point>();
            int iteration = 0;
            int condition = 0;

            // Get Input

            if (!DA.GetData(0, ref mesh)) { return; }
            if (!DA.GetData(1, ref iteration)) { return; }
            if (!DA.GetData(2, ref condition)) { return; }

            // Core of the component
            LaplacianSmoothing.Core_NotWeighted(mesh, iteration, condition, out HeMesh<Euc.Point> otherMesh);

            // Set Output
            DA.SetData(0, otherMesh);

        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_K.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid
        {
            get { return new Guid("{4A536B05-C77B-4E9F-A5E5-D14FBF9590E9}"); }
        }

        /// <inheritdoc cref="GH_K.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <inheritdoc cref="GH_K.GH_DocumentObject.Exposure"/>
        public override GH_K.GH_Exposure Exposure
        {
            get { return GH_K.GH_Exposure.primary; }
        }

        #endregion
    }
}