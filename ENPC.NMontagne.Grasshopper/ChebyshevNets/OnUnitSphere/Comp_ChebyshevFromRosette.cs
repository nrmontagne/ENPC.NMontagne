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
    /// A grasshopper component computing a Chebyshev net on the unity sphere from a rosette condition.
    /// </summary>
    public class Comp_ChebyshevFromRosette : GH_K.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initializes a new instance of the <see cref="Comp_ChebyshevFromRosette"/> class.
        /// </summary>
        public Comp_ChebyshevFromRosette()
          : base("Chebyshev Rosette", "Cheb. R",
              "Computes a Chebyshev net on the unit sphere from a rosette condition.",
              LibrarySettings.Otter.Name, LibrarySettings.Otter.ChebyshevNets.Name)
        {
        }

        #endregion

        #region Override : GH_Component

        /// <inheritdoc cref="GH_K.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_K.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Param_Point(), "Singularity Point", "S", "The singularity point corresponding to a degenerate dual condition.", GH_K.GH_ParamAccess.item);
            pManager.AddParameter(new Param_Point(), "First Ring", "R", "The list of points corresponding to the ring following the singularity.", GH_K.GH_ParamAccess.list);
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
            Point centre = new Point();
            List<Point> ring = new List<Point>();
            int iteration = 0;

            // Get Input

            if (!DA.GetData(0, ref centre)) { return; }
            if (!DA.GetDataList(1, ring)) { return; }
            if (!DA.GetData(2, ref iteration)) { return; }

            // Core of the component
            ChebyshevOnUnitSphere.Core_Rosette(centre, ring, iteration, out HeMesh<Point> mesh);

            // Set Output
            DA.SetData(0, mesh);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_K.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid
        {
            get { return new Guid("{E8150E3F-574F-4483-A837-9B3297457C3D}"); }
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
