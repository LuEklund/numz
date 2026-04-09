const std = @import("std");
const vec = @import("vector.zig");
const Mat4x4 = @import("matrix.zig").@"4x4";
const Quat = @import("quaternion.zig").Hamiltonian;

pub fn Rotor(T: type) type {
    return struct {
        scalar: T,
        xy: T,
        yz: T,
        zx: T,

        pub const identity: @This() = .{ .scalar = 0, .xy = 0, .yz = 0, .zx = 0 };

        pub fn new(scalar: T, xy: T, yz: T, zx: T) @This() {
            return .{ .scalar = scalar, .xy = xy, .yz = yz, .zx = zx };
        }

        pub fn combine(lhs: @This(), rhs: @This()) @This() {
            var result: @This() = undefined;
            result.scalar = lhs.scalar * rhs.scalar - lhs.xy * rhs.xy - lhs.yz * rhs.yz - lhs.zx * rhs.zx;
            result.xy = lhs.scalar * rhs.xy + lhs.xy * rhs.scalar - lhs.yz * rhs.zx + lhs.zx * rhs.yz;
            result.yz = lhs.scalar * rhs.yz + lhs.xy * rhs.zx + lhs.yz * rhs.scalar - lhs.zx * rhs.xy;
            result.zx = lhs.scalar * rhs.zx - lhs.xy * rhs.yz + lhs.yz * rhs.xy + lhs.zx * rhs.scalar;

            return result;
        }

        pub fn reverse(r: @This()) @This() {
            return .{
                .scalar = r.scalar,
                .xy = -r.xy,
                .yz = -r.yz,
                .zx = -r.zx,
            };
        }

        pub fn transform(r: @This(), v: @Vector(3, T)) @Vector(3, T) {
            const S_x: f32 = r.scalar * v[0] + r.xy * v[1] - r.zx * v[2];
            const S_y: f32 = r.scalar * v[1] - r.xy * v[0] + r.yz * v[2];
            const S_z: f32 = r.scalar * v[2] - r.yz * v[1] + r.zx * v[1];
            const S_xyz = r.xy * v[2] + r.zy * v[0] + r.zx * v[1];

            return .{
                S_x * r.scalar + S_y * r.xy + S_xyz * r.yz - S_z * r.zx,
                S_y * r.scalar - S_x * r.xy + S_z * r.yz + S_xyz * r.zx,
                S_z * r.scalar + S_xyz * r.xy - S_y * r.yz + S_x * r.zx,
            };
        }

        pub fn rotate(r: @This()) Mat4x4(T) {
            const new_x: @Vector(3, T) = transform(r, @Vector(3, f32){ 1.0, 0.0, 0.0 });
            const new_y: @Vector(3, T) = transform(r, @Vector(3, f32){ 0.0, 1.0, 0.0 });
            const new_z: @Vector(3, T) = transform(r, @Vector(3, f32){ 0.0, 0.0, 1.0 });

            return .new(.{ new_x[0], new_x[1], new_x[2], 0.0, new_y[0], new_y[1], new_y[2], 0.0, new_z[0], new_z[1], new_z[2], 0.0, 0.0, 0.0, 0.0, 0.0 });
        }

        pub fn fromQuaternion(q: Quat(T)) @This() {
            return .{
                .scalar = q.w,
                .xy = -q.z,
                .yz = -q.x,
                .zx = -q.y,
            };
        }

        pub fn nlerp(lhs: @This(), rhs: @This(), t: f32) @This() {
            const dot: f32 = lhs.scalar * rhs.scalar + lhs.xy * rhs.xy + lhs.yz * rhs.yz + lhs.zx * rhs.zx;
            if (dot < 0.0) {
                rhs.scalar = -rhs.scalar;
                rhs.xy = -rhs.xy;
                rhs.yz = -rhs.yz;
                rhs.zx = -rhs.zx;
            }

            var r: @This() = .{
                .scalar = std.math.lerp(lhs.scalar, rhs.scalar, t),
                .xy = std.math.lerp(lhs.xy, rhs.xy, t),
                .yz = std.math.lerp(lhs.yz, rhs.yz, t),
                .zx = std.math.lerp(lhs.zx, rhs.zx, t),
            };

            const magnitude: f32 = std.math.sqrt(r.scalar * r.scalar + r.xy * r.xy + r.yz * r.yz + r.zx * r.zx);
            r = .{
                .scalar = r.scalar / magnitude,
                .xy = r.xy / magnitude,
                .yz = r.yz / magnitude,
                .zx = r.zx / magnitude,
            };
            return r;
        }

        pub fn slerp(from: @This(), to: @This(), t: f32) @This() {
            var dot: f32 = from.scalar * to.scalar + from.xy * to.xy + from.yz * to.yz + from.zx * to.zx;
            if (dot < 0.0) {
                to.scalar = -to.scalar;
                to.xy = -to.xy;
                to.yz = -to.yz;
                to.zx = -to.zx;
                dot = -dot;
            }

            if (dot > 0.99995) {
                return nlerp(from, to, t);
            }

            const cos_theta = dot;
            const theta: f32 = std.math.acos(cos_theta);
            const from_factor: f32 = std.math.sin((1.0 - t) * theta / std.math.sin(theta));
            const to_factor: f32 = std.math.sin((t * theta) / std.math.sin(theta));

            return .{
                .scalar = from_factor * from.scalar + to_factor * to.scalar,
                .xy = from_factor * from.xy + to_factor * to.xy,
                .yz = from_factor * from.yz + to_factor * to.yz,
                .zx = from_factor * from.zx + to_factor * to.zx,
            };
        }

        pub fn fromVec(from_dir: @Vector(3, T), to_dir: @Vector(3, T)) @This() {
            // This function might be completely incorrect I may need to re review
            from_dir = vec.normalize(from_dir);
            to_dir = vec.normalize(to_dir);

            const halfway: @Vector(3, T) = vec.normalize(from_dir + to_dir);

            const wedge = @Vector(3, T){
                (halfway[0] * from_dir[1]) - (halfway[1] * from_dir[0]),
                (halfway[1] * from_dir[2]) - (halfway[2] * from_dir[1]),
                (halfway[2] * from_dir[0]) - (halfway[0] * from_dir[2]),
            };

            return .{
                .scalar = vec.dot(from_dir, halfway),
                .xy = wedge.x,
                .yz = wedge.y,
                .zx = wedge.z,
            };
        }
    };
}
