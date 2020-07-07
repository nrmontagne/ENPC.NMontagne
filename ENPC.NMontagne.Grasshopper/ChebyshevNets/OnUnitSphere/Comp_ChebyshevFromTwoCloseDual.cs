using System;
using System.Collections.Generic;

using ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;

using GH_K = Grasshopper.Kernel;

using ENPC.McNeel.Grasshopper.Parameters.Euclidean;

using ENPC.NMontagne.Core.CoreFunctions.ChebyshevNets;


namespace ENPC.NMontagne.Grasshopper.VossNets
{
    /// <summary>
    /// A grasshopper component computing a Chebyshev net on the unity sphere from two close dual conditions.
    /// </summary>
    public class Comp_ChebyshevFromTwoCloseDual : GH_K.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initializes a new instance of the <see cref="Comp_ChebyshevFromTwoCloseDual"/> class.
        /// </summary>
        public Comp_ChebyshevFromTwoCloseDual()
          : base("Chebyshev 2 Close Dual", "Cheb. 2CD",
              "Computes a Chebyshev net on the unit sphere from two close dual conditions",
              LibrarySettings.Otter.Name, LibrarySettings.Otter.ChebyshevNets.Name)
        {
        }

        #endregion

        #region Override : GH_Component

        /// <inheritdoc cref="GH_K.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_K.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Param_Point(), "Main Diagonal", "D0", "The list of points corresponding to the first dual condition.", GH_K.GH_ParamAccess.list);
            pManager.AddParameter(new Param_Point(), "Minor Diagonal", "D1", "The list of points corresponding to the second dual condition.", GH_K.GH_ParamAccess.list);
            pManager.AddIntegerParameter("Iteration Count", "I", "The number of iteration to compute ", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_K.GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Param_HeMesh(), "HeMesh", "M", "Mesh representing the Chebyshev net obtained from the input primal conditions.", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.SolveInstance(GH_K.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_K.IGH_DataAccess DA)
        {
            // Declaration, initialization, and instanciation of input variables
            List<Point> d0 = new List<Point>();
            List<Point> d1 = new List<Point>();
            int iteration = 0;

            // Get Input

            if (!DA.GetDataList(0, d0)) { return; }
            if (!DA.GetDataList(1, d1)) { return; }
            if (!DA.GetData(2, ref iteration)) { return; }

            // Core of the component
            ChebyshevOnUnitSphere.Core_FromTwoCloseDual(d0, d1, iteration, out HeMesh<Point> mesh);

            // Set Output
            DA.SetData(0, mesh);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_K.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid
        {
            get { return new Guid("{C6539B95-0F77-462A-B392-5736A2C6A976}"); }
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
