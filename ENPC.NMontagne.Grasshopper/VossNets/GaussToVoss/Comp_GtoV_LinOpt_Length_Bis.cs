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
    /// A grasshopper component computing a Voss net from a Gauss map using a target length for the edges.
    /// </summary>
    public class Comp_GtoV_LinOpt_Length_Bis : GH_K.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initializes a new instance of the <see cref="Comp_GtoV_LinOpt_Length"/> class.
        /// </summary>
        public Comp_GtoV_LinOpt_Length_Bis()
          : base("Gauss To Voss (L2)", "GToV L2",
              "Computes a Voss net from a Gauss map using a target length for the edges.",
              LibrarySettings.Otter.Name, LibrarySettings.Otter.VossNets.Name)
        {
        }

        #endregion

        #region Override : GH_Component

        /// <inheritdoc cref="GH_K.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_K.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Param_Euc.Param_HeMesh(), "Gauss Map", "G", "The Gauss map of the Voss net to compute.", GH_K.GH_ParamAccess.item);
            pManager.AddBooleanParameter("Invert Gauss Curvature ?", "Kp?", "Evaluates whether the Gaussian curvature should be positive or not.", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_K.GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Param_Euc.Param_HeMesh(), "Voss net", "V", "The Voss net generated from the input Gauss map.", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.SolveInstance(GH_K.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_K.IGH_DataAccess DA)
        {
            // Declaration, initialization, and instanciation of input variables
            HeMesh<Point> gaussMap = new HeMesh<Point>();
            bool Kp = false;

            // Get Input

            if (!DA.GetData(0, ref gaussMap)) { return; }
            if (!DA.GetData(1, ref Kp)) { return; }

            // Core of the component
            GaussToVoss_LinOpt.Core_TargetLength_Sparse(gaussMap, Kp, out HeMesh<Point> vossNet);

            // Set Output
            DA.SetData(0, vossNet);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_K.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid
        {
            get { return new Guid("{0A7692D0-C630-4B7E-8FFA-23AACD8E2CB0}"); }
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
