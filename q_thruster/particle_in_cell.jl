using CUDA

# Constants
const pi			= 3.141592653589793			# pi
const eps0		    = 8.85418782E-12			# Permittivity of free space
const mu0			= 1.2566370614E-6			# Vacuum permeability; N-s2/C2
const mu0o4pi		= 1E-7						# mu0/(4*pi)
const inv4piEps0    = 8.9875517873681764e9		# 1/(4*pi*eps0)
const QE			= 1.602176565E-19			# elementary charge (C)
const ME			= 9.10938215E-31			# electron rest mass (kg)

function particleMover(
    px::AbstractFloat, py::AbstractFloat, pz::AbstractFloat,
    vx::AbstractFloat, vy::AbstractFloat, vz::AbstractFloat,
    pq::Int,
    E_realx::AbstractFloat, E_imagx::AbstractFloat, B_realx::AbstractFloat, B_imagx::AbstractFloat,
    E_realy::AbstractFloat, E_imagy::AbstractFloat, B_realy::AbstractFloat, B_imagy::AbstractFloat,
    E_realz::AbstractFloat, E_imagz::AbstractFloat, B_realz::AbstractFloat, B_imagz::AbstractFloat,
    pxBuffer::AbstractFloat, pyBuffer::AbstractFloat, pzBuffer::AbstractFloat,
    vxBuffer::AbstractFloat, vyBuffer::AbstractFloat, vzBuffer::AbstractFloat,
    pqBuffer::Int,
    DestinationsBuffer::Int,
    TransferIndex::Int, KillIndex::Int,
    TransferFlag::Int, KillFlag::Int,
    NumTransfer::Int, NumKill::Int,
    Bxm::Int, Bxp::Int, Bym::Int, Byp::Int, Bzm::Int, Bzp::Int,
    Bxmin::AbstractFloat, Bxmax::AbstractFloat, Bymin::AbstractFloat, Bymax::AbstractFloat, Bzmin::AbstractFloat, Bzmax::AbstractFloat,
    pnum::Int, nB::Int, dt::AbstractFloat, qom::AbstractFloat,
    LX::AbstractFloat, LY::AbstractFloat,
    cosphase::AbstractFloat, sinphase::AbstractFloat,
    L::AbstractFloat, R1::AbstractFloat, R2::AbstractFloat, inv_thresh::AbstractFloat,
    qRPM::AbstractFloat, INTERACTION::Int)
