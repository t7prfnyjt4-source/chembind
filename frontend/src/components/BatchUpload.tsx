import React, { useMemo, useState } from "react";
import { parseSmilesCsv } from "../utils/csvParser";
import { sha256Hex } from "../utils/sha256";

type Props = {
  token: string;
};

export default function BatchUpload({ token }: Props) {
  const [csvText, setCsvText] = useState("");
  const [rows, setRows] = useState<any[]>([]);
  const [jobId, setJobId] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  const preview = useMemo(() => rows.slice(0, 10), [rows]);

  async function onFile(file: File) {
    setError(null);
    setJobId(null);

    const text = await file.text();
    setCsvText(text);

    const parsed = parseSmilesCsv(text);
    setRows(parsed);

    if (parsed.length === 0) {
      setError("No SMILES rows found.");
    }
  }

  async function submit() {
    if (!token || rows.length === 0) return;

    try {
      const key = await sha256Hex(csvText);

      const res = await fetch("/api/batch", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
          "Idempotency-Key": key,
        },
        body: JSON.stringify({
          rows,
          source: "csv",
        }),
      });

      const data = await res.json();
      setJobId(data.job_id);
    } catch (e) {
      setError("Upload failed.");
    }
  }

  return (
    <div>
      <h2>Batch Upload</h2>

      <input
        aria-label="csv-file"
        type="file"
        accept=".csv"
        onChange={(e) => {
          const f = e.target.files?.[0];
          if (f) onFile(f);
        }}
      />

      {rows.length > 0 && (
        <div aria-label="preview">
          <p>
            Preview (first {Math.min(10, rows.length)} of {rows.length})
          </p>
          <table>
            <thead>
              <tr>
                <th>rowId</th>
                <th>smiles</th>
              </tr>
            </thead>
            <tbody>
              {preview.map((r) => (
                <tr key={r.rowId}>
                  <td>{r.rowId}</td>
                  <td>{r.smiles}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      <button
        aria-label="submit"
        onClick={submit}
        disabled={rows.length === 0}
      >
        Submit
      </button>

      {jobId && (
        <p aria-label="job-id">
          job_id: {jobId}
        </p>
      )}

      {error && (
        <p aria-label="error" style={{ color: "red" }}>
          {error}
        </p>
      )}
    </div>
  );
}
