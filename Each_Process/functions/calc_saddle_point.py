from fenics import *
import mshr
import scipy.io
import numpy as np
from matplotlib import pyplot as plt
from scipy.sparse.linalg import eigsh
import scipy.sparse as sp

def calc_saddle_point():
    # Load mesh data from .mat files
    # nodes = scipy.io.loadmat('mesh_nodes.mat')['nodes']
    # elements = scipy.io.loadmat('mesh_elements.mat')['elements']
    # edges = scipy.io.loadmat('mesh_edges.mat')['edges']
    # domain = scipy.io.loadmat('mesh_domain.mat')['domain']

    # # Convert mesh data to FEniCS format
    # mesh = Mesh()
    # editor = MeshEditor()
    # editor.open(mesh, 'triangle', 2, 2)
    # editor.init_vertices(len(nodes))
    # editor.init_cells(len(elements))

    # for i, node in enumerate(nodes):
    #     editor.add_vertex(i, node[0], node[1])

    # for i, element in enumerate(elements):
    #     editor.add_cell(i, int(element[0]) - 1, int(element[1]) - 1, int(element[2]) - 1)

    # editor.close()

    # Define the vertices of the triangle
    x_vertex = 0.5  # Adjust this value if needed
    y_vertex = 1.0  # Adjust this value if needed

    # Create a domain using mshr
    domain = mshr.Polygon([Point(0, 0), Point(1, 0), Point(x_vertex, y_vertex)])
    mesh = mshr.generate_mesh(domain, 32)  # Adjust resolution if needed

    # Define the eigenvalue solver function
    def get_leading_eigenvalues(MatA, MatB, nreq):
        row_, col_, val_ = MatA.mat().getValuesCSR()
        sA = sp.csr_matrix((val_, col_, row_))  # Sparse matrix
        row_, col_, val_ = MatB.mat().getValuesCSR()
        sB = sp.csr_matrix((val_, col_, row_))  # Sparse matrix

        reverse_eigenvalues, reverse_eigenvectors = eigsh(sB, k=nreq, M=sA, which="LM", return_eigenvectors=True, mode="normal")
        eigenvalues = 1. / reverse_eigenvalues
        eigenvalues.sort()
        return eigenvalues, np.flip(reverse_eigenvectors, 1)

    # Define the function to calculate eigenvalue and eigenfunction
    def calc_eigen(mesh, kth, dim_fem):
        def bdry(x, on_boundary): return on_boundary
        V = FunctionSpace(mesh, 'CG', dim_fem)
        bc = DirichletBC(V, Constant(0.0), bdry)

        u = TrialFunction(V)
        v = TestFunction(V)
        a = dot(grad(u), grad(v)) * dx
        b = dot(u, v) * dx

        L = v * dx

        A, _ = assemble_system(a, L, bc)
        B = assemble(b)
        bc.zero(B)

        MatA = as_backend_type(A)
        MatB = as_backend_type(B)

        eigenvalues, eigenvectors = get_leading_eigenvalues(MatA, MatB, kth + 1)

        rx0 = eigenvectors[:, kth - 1]
        u = Function(V)
        u.vector()[:] = rx0
        uu = assemble(dot(u, u) * dx)
        u = u / np.sqrt(uu)
        lam = eigenvalues[kth - 1]

        return u, lam

    dim_fem = 2
    k = 1

    # Calculate eigenfunction u_{i,h}
    u_i_h, eigenvalue = calc_eigen(mesh, k, dim_fem)

    # Step 1: Solve div(p) + u_i_h = 0 for p
    RT = FunctionSpace(mesh, 'RT', 1)
    p = TrialFunction(RT)
    q = TestFunction(RT)

    a1 = dot(p, q) * dx
    L1 = -u_i_h * div(q) * dx

    p_h = Function(RT)
    solve(a1 == L1, p_h)

    # Step 2: Solve (p_h, q_h) + (g_i_h, div(q_h)) = 0 for g_i_h
    V = FunctionSpace(mesh, 'CG', 1)
    g = TrialFunction(V)
    f = TestFunction(V)

    a2 = dot(grad(g), grad(f)) * dx
    L2 = -div(p_h) * f * dx

    g_i_h = Function(V)
    solve(a2 == L2, g_i_h)


    # Plot the solutions
    plt.figure()
    plot(p_h)
    plt.title("Solution p_h")
    plt.show()

    plt.figure()
    plot(g_i_h)
    plt.title("Solution g_i_h")
    plt.show()

    # Convert solutions to numpy arrays
    p_h_matrix = p_h.vector().get_local().reshape((-1, 1))
    g_i_h_matrix = g_i_h.vector().get_local().reshape((-1, 1))

    # Save results to .mat files for MATLAB
    # scipy.io.savemat('p_h_matrix.mat', {'p_h_matrix': p_h_matrix})
    # scipy.io.savemat('g_i_h_matrix.mat', {'g_i_h_matrix': g_i_h_matrix})

# Call the function if the script is executed
if __name__ == "__main__":
    calc_saddle_point()
