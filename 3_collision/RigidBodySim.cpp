#include "RigidBodySim.h"

/////////////////////////////////////
//////// EX 2 - PROBLEM 1 ///////////
/////////////////////////////////////
bool RigidBodySim::advance() {
    // compute the collision detection
    m_collisionDetection.computeCollisionDetection(m_broadPhaseMethod, m_narrowPhaseMethod, m_eps);

    // apply forces (only gravity in this case)
    for (auto &o : m_objects) {
        o.applyForceToCOM(m_gravity);
    }

    for (auto &o : m_objects) {        
        // integrate velocities
        o.setLinearMomentum(o.getLinearMomentum() + m_dt * o.getForce());
        o.setAngularMomentum(o.getAngularMomentum() + m_dt * o.getTorque());
        o.resetForce();
        o.resetTorque();

        // integrate position
        o.setPosition(o.getPosition() + m_dt * o.getLinearVelocity());

        // integrate rotation
        // angular velocity
        Eigen::Vector3d w = o.getAngularVelocity();

        // update orientation
        switch (m_method) {
        case 0: {
            // matrix-based angular velocity
            Eigen::Matrix3d r = o.getRotationMatrix();
            Eigen::Matrix3d W;

            // skew-matrix (row-wise)
            W << 0, -w.z(), w.y(),
                w.z(), 0, -w.x(),
                -w.y(), w.x(), 0;

            r = r + m_dt * W * r;

            // orthogonalize rotation matrix to show issue
            // https://math.stackexchange.com/questions/3292034/normalizing-a-rotation-matrix
            // https://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem
            // idea is to find the nearest orthogonral matrix by SVD
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(r, Eigen::ComputeFullU | Eigen::ComputeFullV);
            r = svd.matrixU() * svd.matrixV().transpose();

            o.setRotation(r);
            break;
        }
        default: {
            // quaternion-based
            Eigen::Quaterniond wq;
            wq.w() = 0;
            wq.vec() = w;

            Eigen::Quaterniond q = o.getRotation();
            Eigen::Quaterniond dq = wq * q;
            Eigen::Quaterniond new_q;
            new_q.w() = q.w() + 0.5 * m_dt * dq.w();
            new_q.vec() = q.vec() + 0.5 * m_dt * dq.vec();
            o.setRotation(new_q.normalized());
            break;
        }
        }
    }

    // advance time
    m_time += m_dt;
    m_step++;

    return false;
}