using System;
using System.Collections.Generic;

using ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;

using GH_K = Grasshopper.Kernel;

using ENPC.McNeel.Grasshopper.Parameters.Euclidean;

using ENPC.NMontagne.Core.CoreFunctions.ChebyshevNets;


namespace ENPC.NMontagne.Grasshopper.ChebyshevNet.OnUnitSphere
{
    /// <summary>
    /// A grasshopper component omputing a Chebyshev net on the unity sphere from two primal conditions.
    /// </summary>
    public class Comp_ChebyshevFromTwoPrimal : GH_K.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initializes a new instance of the <see cref="Comp_ChebyshevFromTwoPrimal"/> class.
        /// </summary>
        public Comp_ChebyshevFromTwoPrimal()
          : base("Chebyshev 2 Primal", "Cheb. 2P",
              "Computes a Chebyshev net on the unit sphere from two primal conditions.",
              LibrarySettings.Otter.Name, LibrarySettings.Otter.ChebyshevNets.Name)
        {
        }

        #endregion

        #region Override : GH_Component

        /// <inheritdoc cref="GH_K.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_K.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Param_Point(), "U Normals", "U", "The list of points corresponding to the primal condition in the first direction U.", GH_K.GH_ParamAccess.list);
            pManager.AddParameter(new Param_Point(), "V Normals", "V", "The list of points corresponding to the primal condition in the first direction V.", GH_K.GH_ParamAccess.list);
        }

        /// <inheritdoc cref="GH_K.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_K.GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Param_HeMesh(), "HeMesh", "M", "The mesh representing the Chebyshev net obtained from the input primal conditions.", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.SolveInstance(GH_K.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_K.IGH_DataAccess DA)
        {
            // Declaration, initialization, and instanciation of input variables
            List<Point> u = new List<Point>();
            List<Point> v = new List<Point>();

            // Get Input

            if (!DA.GetDataList(0, u)) { return; }
            if (!DA.GetDataList(1, v)) { return; }

            // Core of the component
            ChebyshevOnUnitSphere.Core_FromTwoPrimal(u, v, out HeMesh<Point> mesh);

            // Set Output
            DA.SetData(0, mesh);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_K.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid
        {
            get { return new Guid("{E2016358-6EA4-423B-ABCC-68910EBAA3F7}"); }
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
