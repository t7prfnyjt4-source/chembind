import { useState } from "react";
import { apiFetch } from "../api/client";
import { useDockingJob } from "../hooks/useDockingJob";
import DockingViewer from "./DockingViewer";

export default function DockingScreen() {
  const [smiles, setSmiles] = useState("");
  const [pdbFile, setPdbFile] = useState<File | null>(null);
  const [centerX, setCenterX] = useState("0");
  const [centerY, setCenterY] = useState("0");
  const [centerZ, setCenterZ] = useState("0");
  const [exhaustiveness, setExhaustiveness] = useState(8);
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [jobId, setJobId] = useState<string | null>(null);

  const job = useDockingJob(jobId);

  const handleSubmit = async () => {
    if (!smiles.trim() || !pdbFile) return;
    setSubmitting(true);
    setError(null);

    try {
      const formData = new FormData();
      formData.append("smiles", smiles.trim());
      formData.append("pdb_file", pdbFile);
      formData.append("center_x", centerX);
      formData.append("center_y", centerY);
      formData.append("center_z", centerZ);
      formData.append("exhaustiveness", String(exhaustiveness));

      const resp = await apiFetch<{ jobId: string }>("/api/docking/jobs", {
        method: "POST",
        body: formData,
        headers: {},
      });
      setJobId(resp.jobId);
    } catch (e: any) {
      setError(e?.message || "Submission failed");
    } finally {
      setSubmitting(false);
    }
  };

  const statusBadge = (status: string) => {
    const colors: Record<string, string> = {
      queued: "#888",
      preparing: "#e8a317",
      docking: "#3b82f6",
      completed: "#22c55e",
      failed: "#ef4444",
    };
    return (
      <span style={{
        padding: "2px 8px",
        borderRadius: 4,
        background: colors[status] || "#888",
        color: "#fff",
        fontSize: 12,
        fontWeight: 600,
      }}>
        {status}
      </span>
    );
  };

  return (
    <div>
      <div style={{ fontWeight: 600, marginBottom: 8 }}>Molecular Docking</div>

      {!jobId && (
        <div style={{ display: "flex", flexDirection: "column", gap: 8, maxWidth: 400 }}>
          <input
            type="text"
            value={smiles}
            onChange={(e) => setSmiles(e.target.value)}
            placeholder="Ligand SMILES (e.g. c1ccccc1)"
            style={{ padding: "6px 10px", fontFamily: "monospace" }}
          />
          <div>
            <label style={{ fontSize: 13, display: "block", marginBottom: 4 }}>PDB File:</label>
            <input
              type="file"
              accept=".pdb"
              onChange={(e) => setPdbFile(e.target.files?.[0] || null)}
            />
          </div>
          <div style={{ display: "flex", gap: 8 }}>
            <div>
              <label style={{ fontSize: 12 }}>Center X</label>
              <input type="number" value={centerX} onChange={(e) => setCenterX(e.target.value)} style={{ width: 60 }} />
            </div>
            <div>
              <label style={{ fontSize: 12 }}>Center Y</label>
              <input type="number" value={centerY} onChange={(e) => setCenterY(e.target.value)} style={{ width: 60 }} />
            </div>
            <div>
              <label style={{ fontSize: 12 }}>Center Z</label>
              <input type="number" value={centerZ} onChange={(e) => setCenterZ(e.target.value)} style={{ width: 60 }} />
            </div>
          </div>
          <div>
            <label style={{ fontSize: 12 }}>Exhaustiveness: {exhaustiveness}</label>
            <input
              type="range"
              min={1}
              max={32}
              value={exhaustiveness}
              onChange={(e) => setExhaustiveness(Number(e.target.value))}
              style={{ width: "100%" }}
            />
          </div>
          <button onClick={handleSubmit} disabled={submitting || !smiles.trim() || !pdbFile}>
            {submitting ? "Submitting…" : "Start Docking"}
          </button>
        </div>
      )}

      {error && <div style={{ color: "#c00", marginTop: 8, fontSize: 14 }}>{error}</div>}

      {job && (
        <div style={{ marginTop: 12 }}>
          <div style={{ display: "flex", gap: 8, alignItems: "center", marginBottom: 8 }}>
            <span style={{ fontSize: 13 }}>Job: {jobId}</span>
            {statusBadge(job.status)}
          </div>

          {job.status === "completed" && (
            <DockingViewer jobId={jobId!} poseCount={job.poseCount || 0} bestScore={job.bestScore} />
          )}

          {job.status === "failed" && (
            <div style={{ color: "#c00", fontSize: 14 }}>Error: {job.error}</div>
          )}
        </div>
      )}
    </div>
  );
}
