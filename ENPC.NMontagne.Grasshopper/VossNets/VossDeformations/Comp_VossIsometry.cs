using System;
using System.Collections.Generic;

using ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;

using GH_K = Grasshopper.Kernel;

using Param_Euc = ENPC.McNeel.Grasshopper.Parameters.Euclidean;

using ENPC.NMontagne.Core.CoreFunctions.VossNets;


namespace ENPC.NMontagne.Grasshopper.VossNets
{
    /// <summary>
    /// A grasshopper component computing the isometric deformation of a Voss net.
    /// </summary>
    public class Comp_VossIsometry : GH_K.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initializes a new instance of the <see cref="Comp_VossIsometry"/> class.
        /// </summary>
        public Comp_VossIsometry()
          : base("Voss Isometry", "Voss Iso.",
              "Computes the isometric deformation of a Voss net.",
              LibrarySettings.Otter.Name, LibrarySettings.Otter.VossNets.Name)
        {
        }

        #endregion

        #region Override : GH_Component

        /// <inheritdoc cref="GH_K.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_K.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Param_Euc.Param_HeMesh(), "Voss Net", "M", "The Voss net to deform.", GH_K.GH_ParamAccess.item);
            pManager.AddNumberParameter("Factor", "F", "The paramater for the isometry of Voss nets.", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_K.GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Param_Euc.Param_HeMesh(), "Voss net", "V", "The deformet Voss net", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.SolveInstance(GH_K.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_K.IGH_DataAccess DA)
        {
            // Declaration, initialization, and instanciation of input variables
            HeMesh<Point> gaussMap = new HeMesh<Point>();
            double factor = 0.0;

            // Get Input

            if (!DA.GetData(0, ref gaussMap)) { return; }
            if (!DA.GetData(1, ref factor)) { return; }

            // Core of the component
            VossDeformations.Core_Isometric(gaussMap, factor, out HeMesh<Point> vossNet);

            // Set Output
            DA.SetData(0, vossNet);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_K.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid
        {
            get { return new Guid("{921A3563-7D99-4C3C-AEB5-7F1DF3F7A8E2}"); }
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
