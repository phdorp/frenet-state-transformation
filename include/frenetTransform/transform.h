#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <Eigen/Core>
#include <memory>

#include "frenetTransform/internal/math.h"
#include "frenetTransform/path.h"
#include "frenetTransform/points.h"

namespace FrenetTransform {
/**
 * @brief Transformation between Cartesian and Frenet frame.
 * Use the Path properties to implement the transformation independent from the
 * Path implementation.
 *
 * @tparam NumQueries number of query points with -1 for dynamic point number.
 */
template <int NumQueries = Eigen::Dynamic> class Transform {
public:
  using ArrayQueries = Eigen::Array<double, NumQueries, 1>;

  Transform() = default;

  /**
   * @brief Construct a new Transform object from a given Path.
   *
   * @param path defines the Frenet frame.
   */
  Transform(const std::shared_ptr<Path<NumQueries>> path) : m_path{path} {}

  /**
   * @brief Transform Cartesian positions to Frenet positions.
   * Projects the query points onto the path.
   * Dertermines the singed lengths along the path from the path origin to the
   * projections. Determines the signed shortest distances to the query point.
   *
   * @param posCartes query points in Cartesian coordinates.
   * @return Points<NumQueries> resut points in Frenet coordinates.
   */
  Points<NumQueries> posFrenet(const Points<NumQueries> &posCartes) const {
    // lengths from path origin to Cartesian positions
    const auto lengths{m_path->lengths(posCartes)};
    // next points on path to Cartesian positions
    const auto posPath{m_path->operator()(lengths)};
    // vectors from query points to path points
    const auto posDiff{posCartes - posPath};
    // normal vectors on path
    const auto normals{m_path->normal(lengths)};
    return {lengths, normals * posDiff};
  }

  /**
   * @brief Transform Frenet positions to Cartesian positions.
   *
   * @param posFrenet query points in Frenet coordinates.
   * @return Points<NumQueries> result points in Cartesian coordinates.
   */
  Points<NumQueries> posCartes(const Points<NumQueries> &posFrenet) const {
    // positions along the path at the signed lenghts from the origin
    const auto posPath{m_path->operator()(posFrenet.x())};
    // normals along the path pointing toward or away from the Cartesian point
    const auto normals{m_path->normal(posFrenet.x())};
    return posPath + normals * posFrenet.y();
  }

  /**
   * @brief Transform Cartesian velocitiese to Frenet velocities.
   *
   * @param velCartes query velocities in Cartesian coordinates.
   * @param posFrenet positions corresponding to velocities in Frenet
   * coordinates.
   * @return Points<NumQueries> result velocities in Frenet coordinates.
   */
  Points<NumQueries> velFrenet(const Points<NumQueries> &velCartes,
                               const Points<NumQueries> &posFrenet) const {
    // transformation matrices from Cartesian to Frenet frames at given Frenet
    // frame positions
    const auto velTransformsInv{transformInv(velTransform(posFrenet))};
    // matrix-vector multiplications of transformations with velocities
    return {velTransformsInv(0, 0) * velCartes.x() +
                velTransformsInv(0, 1) * velCartes.y(),
            velTransformsInv(1, 0) * velCartes.x() +
                velTransformsInv(1, 1) * velCartes.y()};
  }

  /**
   * @brief Tranform Frenet velocities to Cartesian velocities.
   *
   * @param velFrenet query velocities in Frenet coordinates.
   * @param posFrenet positions corresponding to velocities in Frenet
   * coordinates.
   * @return Points<NumQueries> result velocities in Cartesian coordinates.
   */
  Points<NumQueries> velCartes(const Points<NumQueries> &velFrenet,
                               const Points<NumQueries> &posFrenet) const {
    // transformation matrices from Frenet to Cartesian frames at given Frenet
    // frame positions
    const auto velTransforms{velTransform(posFrenet)};
    // matrix-vector multiplications of transformations with velocities
    return {velTransforms(0, 0) * velFrenet.x() +
                velTransforms(0, 1) * velFrenet.y(),
            velTransforms(1, 0) * velFrenet.x() +
                velTransforms(1, 1) * velFrenet.y()};
  }

  /**
   * @brief Transform Cartesian accelerations to Frenet accelerations.
   *
   * @param accCartes query accelerations in Cartesian coordinates.
   * @param velFrenet velocities corresponding to accelerations in Frenet
   * coordinates.
   * @param posFrenet positions corresponding to accelerations in Frenet
   * coordinates.
   * @return Points<NumQueries> result accelerations in Frenet coordinates.
   */
  Points<NumQueries> accFrenet(const Points<NumQueries> &accCartes,
                               const Points<NumQueries> &velFrenet,
                               const Points<NumQueries> &posFrenet) const {
    // transformation matrices to determine Cartesian acceleration induced by
    // Frenet velocities
    const auto accTransforms{accTransform(velFrenet, posFrenet)};
    // difference between query accelerations and velocity-induced accelerations
    const Points<NumQueries> accDiff{
        accCartes.x() - accTransforms(0, 0) * velFrenet.x() -
            accTransforms(0, 1) * velFrenet.y(),
        accCartes.y() - accTransforms(1, 0) * velFrenet.x() -
            accTransforms(1, 1) * velFrenet.y()};
    // transformation matrices from Cartesian to Frenet frames at given Frenet
    // frame positions
    const auto velTransformsInv{transformInv(velTransform(posFrenet))};
    // matrix-vector multiplications of transformations with velocity
    // differences
    return {velTransformsInv(0, 0) * accDiff.x() +
                velTransformsInv(0, 1) * accDiff.y(),
            velTransformsInv(1, 0) * accDiff.x() +
                velTransformsInv(1, 1) * accDiff.y()};
  }

  /**
   * @brief Transform Frenet accelerations to Cartesian accelerations.
   *
   * @param accFrenet query accelerations in Frenet coordinates.
   * @param velFrenet velocities corresponding to accelerations in Frenet
   * coordinates.
   * @param posFrenet positions corresponding to accelerations in Frenet
   * coordinates.
   * @return Points<NumQueries> result accelerations in Cartesian coordinates.
   */
  Points<NumQueries> accCartes(const Points<NumQueries> &accFrenet,
                               const Points<NumQueries> &velFrenet,
                               const Points<NumQueries> &posFrenet) const {
    // transformation matrices to determine Cartesian acceleration induced by
    // Frenet velocities
    const auto accTransforms{accTransform(velFrenet, posFrenet)};
    // transformation matrices from Frenet to Cartesian frames at given Frenet
    // frame positions
    const auto velTransforms{velTransform(posFrenet)};
    // sum of velocity-induced accelerations and transformed Frenet
    // accelerations
    return {velTransforms(0, 0) * accFrenet.x() +
                velTransforms(0, 1) * accFrenet.y() +
                accTransforms(0, 0) * velFrenet.x() +
                accTransforms(0, 1) * velFrenet.y(),
            velTransforms(1, 0) * accFrenet.x() +
                velTransforms(1, 1) * accFrenet.y() +
                accTransforms(1, 0) * velFrenet.x() +
                accTransforms(1, 1) * velFrenet.y()};
  }

protected:
  std::shared_ptr<Path<NumQueries>> m_path; /**< Store path. */

private:
  /**
   * @brief Transformation matrices from Frenet to Cartesian frame.
   *
   * @param posFrenet query positions in Frenet frame.
   * @return Eigen::Array<ArrayQueries, 2, 2> matrices at "posFrenet".
   */
  Eigen::Array<ArrayQueries, 2, 2>
  velTransform(const Points<NumQueries> &posFrenet) const {
    const auto tangents{m_path->tangent(posFrenet.x())};
    const auto normals{m_path->normal(posFrenet.x())};
    const auto curvs{m_path->angle1(posFrenet.x())};
    return {{tangents.x() * (1 - curvs * posFrenet.y()), normals.x()},
            {tangents.y() * (1 - curvs * posFrenet.y()), normals.y()}};
  }

  /**
   * @brief Transformation matrices from Frenet velocities to Cartesian
   * velocity-induced accelerations.
   *
   * @param velFrenet velocities in Frenet frame.
   * @param posFrenet positions in Frenet frame.
   * @return Eigen::Array<ArrayQueries, 2, 2> transformations at "velFrenet" and
   * "posFrenet".
   */
  Eigen::Array<ArrayQueries, 2, 2>
  accTransform(const Points<NumQueries> &velFrenet,
               const Points<NumQueries> &posFrenet) const {
    const auto tangents{m_path->tangent(posFrenet.x())};
    const auto normals{m_path->normal(posFrenet.x())};
    // scaling factor due to lateral position from path and path curvature
    const auto curvs{m_path->angle1(posFrenet.x())};
    const auto latScale{1 - curvs * posFrenet.y()};
    // scaling factor derivatvie
    const auto curv1s{m_path->angle2(posFrenet.x())};
    const auto latScaleDer{curv1s * velFrenet.x() * posFrenet.y() +
                           curvs * velFrenet.y()};
    return {{normals.x() * curvs * latScale * velFrenet.x() -
                 tangents.x() * latScaleDer,
             -curvs * tangents.x() * velFrenet.x()},
            {normals.y() * curvs * latScale * velFrenet.x() -
                 tangents.y() * latScaleDer,
             -curvs * tangents.y() * velFrenet.x()}};
  }
};
}; // namespace FrenetTransform

#endif