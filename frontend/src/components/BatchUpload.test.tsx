import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import BatchUpload from "./BatchUpload";

/**
 * Proper crypto.subtle.digest mock (JSDOM safe)
 * We spy on the existing subtle.digest instead of overwriting crypto.
 */
function mockCryptoDigest(hex64: string) {
  const bytes = new Uint8Array(
    hex64.match(/.{1,2}/g)!.map((b) => parseInt(b, 16))
  );
  const buf = bytes.buffer;

  vi.spyOn(global.crypto.subtle, "digest").mockResolvedValue(buf);
}

describe("BatchUpload", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("loads CSV, shows preview, submits, displays job_id", async () => {
    mockCryptoDigest(
      "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"
    );

    const fetchMock = vi.fn(async () => {
      return new Response(JSON.stringify({ job_id: "job_123" }), {
        status: 200,
        headers: { "Content-Type": "application/json" },
      });
    });

    // @ts-expect-error - test overrides fetch with a mock
    global.fetch = fetchMock;

    render(<BatchUpload token="TOKEN_ABC" />);

    const csv = `smiles
CCO
CCC
`;

    const file = new File([csv], "x.csv", { type: "text/csv" });

    const input = screen.getByLabelText("csv-file");
    fireEvent.change(input, { target: { files: [file] } });

    // Wait for preview to render
    await screen.findByLabelText("preview");

    expect(screen.getByText(/first 2 of 2/i)).toBeInTheDocument();
    expect(screen.getByText("CCO")).toBeInTheDocument();
    expect(screen.getByText("CCC")).toBeInTheDocument();

    fireEvent.click(screen.getByLabelText("submit"));

    await waitFor(() => {
      expect(fetchMock).toHaveBeenCalledTimes(1);
    });

    const [, opts] = fetchMock.mock.calls[0];

    expect(opts.method).toBe("POST");
    expect(opts.headers.Authorization).toBe("Bearer TOKEN_ABC");
    expect(opts.headers["Idempotency-Key"]).toMatch(/^[0-9a-f]{64}$/);
    expect(opts.headers["Content-Type"]).toBe("application/json");

    const body = JSON.parse(opts.body);
    expect(body.source).toBe("csv");
    expect(body.rows.length).toBe(2);
    expect(body.rows[0].smiles).toBe("CCO");

    const jobEl = await screen.findByLabelText("job-id");
    expect(jobEl).toHaveTextContent("job_123");
  });
});
