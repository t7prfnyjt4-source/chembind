import { useState, useEffect } from "react";
import { apiFetch } from "../api/client";

interface PoseAtom {
  name: string;
  x: number;
  y: number;
  z: number;
  element: string;
}

interface Pose {
  pose_id: number;
  score: number;
  atoms: PoseAtom[];
  interactions?: { residue: string; type: string }[];
}

interface PosesResponse {
  poses: Pose[];
}

const ELEMENT_COLORS: Record<string, string> = {
  C: "#333", N: "#3050F8", O: "#FF0D0D", S: "#FFFF30",
  P: "#FF8000", H: "#ccc",
};

function projectAtoms(atoms: PoseAtom[], w: number, h: number) {
  if (!atoms.length) return [];
  const xs = atoms.map(a => a.x), ys = atoms.map(a => a.y);
  const minX = Math.min(...xs), maxX = Math.max(...xs);
  const minY = Math.min(...ys), maxY = Math.max(...ys);
  const rX = maxX - minX || 1, rY = maxY - minY || 1;
  const scale = Math.min((w - 40) / rX, (h - 40) / rY);
  const cx = w / 2, cy = h / 2, mx = (minX + maxX) / 2, my = (minY + maxY) / 2;
  return atoms.map(a => ({
    sx: cx + (a.x - mx) * scale,
    sy: cy - (a.y - my) * scale,
    element: a.element,
    color: ELEMENT_COLORS[a.element] || "#888",
  }));
}

export default function DockingViewer({
  jobId,
  poseCount,
  bestScore,
}: {
  jobId: string;
  poseCount: number;
  bestScore?: number;
}) {
  const [poses, setPoses] = useState<Pose[]>([]);
  const [selectedPose, setSelectedPose] = useState(0);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    apiFetch<PosesResponse>(`/api/docking/jobs/${jobId}/poses`)
      .then((resp) => {
        setPoses(resp.poses);
        setLoading(false);
      })
      .catch(() => setLoading(false));
  }, [jobId]);

  if (loading) return <div style={{ opacity: 0.7 }}>Loading poses…</div>;
  if (!poses.length) return <div>No poses found.</div>;

  const current = poses[selectedPose];
  const projected = projectAtoms(current.atoms, 300, 300);

  return (
    <div>
      <div style={{ fontSize: 13, marginBottom: 8, padding: "6px 10px", background: "#fff3cd", borderRadius: 4 }}>
        For research purposes only. Not for clinical use.
      </div>

      <table style={{ width: "100%", borderCollapse: "collapse", fontSize: 13, marginBottom: 12 }}>
        <thead>
          <tr style={{ borderBottom: "1px solid #555" }}>
            <th style={{ textAlign: "left", padding: 4 }}>Pose</th>
            <th style={{ textAlign: "left", padding: 4 }}>Score (kcal/mol)</th>
            <th style={{ padding: 4 }}></th>
          </tr>
        </thead>
        <tbody>
          {poses.map((p, i) => (
            <tr key={p.pose_id} style={{ borderBottom: "1px solid #eee", background: i === selectedPose ? "#e8f0fe" : "transparent" }}>
              <td style={{ padding: 4 }}>#{p.pose_id}</td>
              <td style={{ padding: 4 }}>{p.score.toFixed(2)}</td>
              <td style={{ padding: 4 }}>
                <button onClick={() => setSelectedPose(i)} style={{ fontSize: 12 }}>View</button>
              </td>
            </tr>
          ))}
        </tbody>
      </table>

      <div style={{ fontWeight: 600, fontSize: 13, marginBottom: 4 }}>
        Pose #{current.pose_id} — {current.score.toFixed(2)} kcal/mol
      </div>

      <svg width={300} height={300} style={{ border: "1px solid #555", borderRadius: 4, background: "#fafafa" }}>
        {projected.map((p, i) => (
          <circle key={i} cx={p.sx} cy={p.sy} r={5} fill={p.color} />
        ))}
      </svg>

      {current.interactions && current.interactions.length > 0 && (
        <div style={{ marginTop: 8, fontSize: 13 }}>
          <div style={{ fontWeight: 600, marginBottom: 4 }}>Interactions:</div>
          {current.interactions.map((int, i) => (
            <span key={i} style={{
              display: "inline-block", padding: "2px 6px", margin: 2,
              borderRadius: 4, background: "#e0e7ff", fontSize: 11,
            }}>
              {int.residue}: {int.type}
            </span>
          ))}
        </div>
      )}

      {(!current.interactions || current.interactions.length === 0) && (
        <div style={{ marginTop: 4, fontSize: 12, opacity: 0.6 }}>No interactions detected</div>
      )}
    </div>
  );
}
